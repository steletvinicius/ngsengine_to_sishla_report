import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import time
import zipfile
from io import BytesIO

# Importando as funções dos seus módulos locais
from xml_read_engine_functions import *
from format_dfs_to_import import *
from readenginetar import *
from pyard_functions import *
from mac_request_api import *

# Criação de Abas na Interface
tab_process, tab_dash = st.tabs(["🚀 Processamento", "📊 Dashboard de Qualidade"])

with tab_process:

    st.set_page_config(page_title="Gerador de Relatórios SISHLA/REDOME a partir de Arquivos XML do NGSEngine", layout="wide")

    st.title("Ferramenta de Importação NGS - NMDP & LIMS")
    st.markdown("Faça o upload dos arquivos XML gerados pelo sequenciador para processar as tipagens e gerar os relatórios do SISHLA e REDOME.")

    # Componentes de Upload atualizados para aceitar XML e ZIP
    uploaded_files = st.file_uploader("Selecione o arquivo export TAR do NGSEngine (XML ou ZIP)", type=['xml', 'zip'], accept_multiple_files=True)
    uploaded_dmr = st.file_uploader("Opcional: Arquivo IL_DMR (CSV) para o REDOME", type=['csv'])

    # Botão para iniciar o processamento
    if st.button("Processar Arquivos") and uploaded_files:
        
        # --- NOVO BLOCO: FILA DE ARQUIVOS E DESCOMPACTAÇÃO EM MEMÓRIA ---
        xml_files_to_process = []
        
        for uploaded_file in uploaded_files:
            if uploaded_file.name.endswith('.zip'):
                # Abre o zip diretamente da memória
                with zipfile.ZipFile(uploaded_file, 'r') as z:
                    for filename in z.namelist():
                        # Filtra apenas os XMLs e ignora arquivos ocultos de sistema (ex: MacOS)
                        if filename.endswith('.xml') and not filename.startswith('__MACOSX'):
                            file_data = z.open(filename) # Retorna um objeto file-like compatível com o extract_metrics_xml
                            xml_files_to_process.append({
                                'name': filename.split("/")[-1], # Pega apenas o nome do arquivo, ignorando as pastas internas do zip
                                'file_object': file_data
                            })
            elif uploaded_file.name.endswith('.xml'):
                xml_files_to_process.append({
                    'name': uploaded_file.name,
                    'file_object': uploaded_file
                })
                
        if not xml_files_to_process:
            st.warning("Nenhum arquivo XML válido foi encontrado dentro dos envios.")
            st.stop()
        # ----------------------------------------------------------------

        # Barra de progresso para a interface
        progress_bar = st.progress(0)
        status_text = st.empty()
        
        # O loop agora itera sobre a fila unificada que criamos
        for i, item in enumerate(xml_files_to_process):
            batch_name = item['name'].split(".")[0]
            status_text.text(f"Processando lote: {batch_name}...")
            
            try:
                # 1. Extração do XML (passando o objeto em memória extraído)
                df = extract_metrics_xml(item['file_object'])
                
                # --- CAPTURA DINÂMICA DA VERSÃO IMGT ---
                imgt_version = '3.62.0'  # Versão padrão de fallback caso algo falhe
                
                if 'imgt_version' in df.columns and not df['imgt_version'].empty:
                    # Obtém o valor bruto da primeira linha (ex: "IMGT/HLA 3.62.0")
                    raw_version = df['imgt_version'].dropna().iloc[0]
                    
                    # Expressão regular para capturar apenas o padrão de números e pontos (ex: 3.62.0)
                    import re
                    match = re.search(r'(\d+\.\d+\.\d+)', str(raw_version))
                    if match:
                        imgt_version = match.group(1)
                
                # Exibe na interface web qual versão foi detectada e carregada
                st.info(f"Lote `{batch_name}`: Base de dados IMGT {imgt_version} detectada e carregada automaticamente.")
                
                # Inicializa o objeto ARD com a versão correta extraída do próprio arquivo
                # Buscará diretamente a pasta correspondente dentro de './py_ard_db'
                ard = init_pyard(imgt_version)
                # ----------------------------------------

                df_typing = df.copy()
                df_quality_metrics = df_typing.copy()

                # Métricas de Qualidade Adicionais
                qm_cols = ['fastq_filename', 'sample_name', 'software_name', 'software_version', "imgt_version", "locus_name",
                        'nmdp_typing_allele1', 'nmdp_typing_allele2', "typing_pgroup", "locus_review_status",
                        "locus_review_level", "locus_first_reviewer_name", "locus_first_reviewer_datetime",
                        "locus_first_reviewer_action", "locus_second_reviewer_name", "locus_second_reviewer_action",
                        "locus_map_percentage", "exon_qm_noise_delta_ston", 'core_qm_noise_delta_ston', 'exon_qm_noise_max',
                        'core_qm_noise_max', "exon_qm_read_depth_min", 'exon_qm_read_depth_max', 'exon_qm_secondbase_min',
                        'core_qm_read_depth_min', 'core_qm_secondbase_min']
                
                df_quality_metrics = df_quality_metrics[qm_cols]

                # 2. Limpeza e Filtros
                df_typing = df_typing.loc[(df_typing['locus_review_status'] == 'Approved') & (~df_typing['sample_name'].str.contains('1044736'))]
                df_typing.drop_duplicates(subset=['sample_name','locus_name'], inplace=True)

                classical_loci = ["HLA-A","HLA-B","HLA-C", "HLA-DRB1", "HLA-DRB3", "HLA-DRB4", "HLA-DRB5",
                                "HLA-DQB1","HLA-DPB1", "HLA-DQA1","HLA-DPA1"]
                df_typing = df_typing.loc[df_typing["locus_name"].isin(classical_loci)]

                # --- NOVO BLOCO: CACHE MULTI-VERSÃO IMGT ---
                import re
                
                # 1. Limpar a coluna imgt_version para extrair apenas o número (ex: "IMGT/HLA 3.62.0" -> "3.62.0")
                def extract_version(text):
                    if pd.notnull(text):
                        match = re.search(r'(\d+\.\d+\.\d+)', str(text))
                        return match.group(1) if match else '3.62.0'
                    return '3.62.0' # Fallback de segurança

                df_typing['imgt_version_clean'] = df_typing['imgt_version'].apply(extract_version)

                # 2. Descobrir quais versões únicas existem neste dataframe
                versoes_unicas = df_typing['imgt_version_clean'].unique()
                st.info(f"Lote `{batch_name}`: Versões IMGT detectadas: {', '.join(versoes_unicas)}")

                # 3. Inicializar o pyard UMA VEZ para cada versão encontrada (Cache)
                ard_cache = {}
                for v in versoes_unicas:
                    ard_cache[v] = init_pyard(v)

                # --- REDUÇÃO PARA 3 CAMPOS USANDO A VERSÃO CORRETA DA LINHA ---
                def aplicar_reducao(row, col_name):
                    valor = row[col_name]
                    versao = row['imgt_version_clean']
                    ard_instancia = ard_cache[versao] # Pega o pyard correto do cache
                    
                    if (valor is not None) and (len(str(valor).split(':')) > 3):
                        return typing_3_fields(ard_instancia, valor)
                    return valor

                # Usamos apply(axis=1) para ter acesso à linha inteira (e poder ler a versão correspondente)
                df_typing['nmdp_typing_allele1'] = df_typing.apply(lambda row: aplicar_reducao(row, 'nmdp_typing_allele1'), axis=1)
                df_typing['nmdp_typing_allele2'] = df_typing.apply(lambda row: aplicar_reducao(row, 'nmdp_typing_allele2'), axis=1)
                
                df_typing['_alelo01+alelo02'] = df_typing["nmdp_typing_allele1"] + "+" + df_typing["nmdp_typing_allele2"]
                
                # --- REQUISIÇÃO MAC API COM A VERSÃO CORRETA DA LINHA ---
                df_typings_to_mac = df_typing.loc[df_typing['_alelo01+alelo02'].str.contains(r"\?", na=False)].copy()
                
                if not df_typings_to_mac.empty:
                    status_text.text(f"Consultando API NMDP MAC para {len(df_typings_to_mac)} alelos...")
                    # Passa o typing_result e a imgt_version_clean específicos daquela linha para a API
                    df_typings_to_mac['_alelo01+alelo02'] = df_typings_to_mac.apply(
                        lambda row: encode_mac(row['typing_result'], row['imgt_version_clean']), axis=1
                    )

                df_typing_complete = df_typing.loc[~df_typing['_alelo01+alelo02'].str.contains(r"\?", na=False)]
                df_typing_complete = pd.concat([df_typing_complete, df_typings_to_mac])

                # 5. Separação de Alelos
                df_typing_complete['_alelo01'] = df_typing_complete['_alelo01+alelo02'].apply(lambda x: str(x).split('+')[0] if x != None else "")
                df_typing_complete['_alelo02'] = df_typing_complete['_alelo01+alelo02'].apply(lambda x: str(x).split("+")[1] if x != None and len(x.split("+"))>1 else "")

                # 6. Formatação SISHLA
                df_sishla = format_df(df_typing_complete, 'sishla')
                df_sishla = df_sishla.rename(columns={"sample_name" : "00_sample_name"})
                df_sishla = df_sishla.reindex(sorted(df_sishla.columns), axis=1)

                # 7. Formatação REDOME
                df_redome = format_df(df_typing_complete, 'redome')
                
                # Se o usuário fez upload da planilha DMR, fazemos o merge
                if uploaded_dmr is not None:
                    df_ils_dmrs = pd.read_csv(uploaded_dmr, dtype=str, encoding="latin", sep=";")
                    df_redome = pd.merge(df_redome, df_ils_dmrs[["Patient", "DMR"]], left_on='sample_name', right_on="Patient", how='left')
                    df_redome = df_redome.rename(columns={"DMR" : "00_DMR"})
                    df_redome = df_redome.reindex(sorted(df_redome.columns), axis=1)

                # --- GERAÇÃO DOS ARQUIVOS PARA DOWNLOAD ---
                st.success(f"Lote {batch_name} processado com sucesso!")
                
                # Função auxiliar para converter DF para CSV em memória
                @st.cache_data
                def convert_df(df):
                    return df.to_csv(sep=';', index=False, encoding='utf-8').encode('utf-8')

                # Organizando os botões de download lado a lado
                # --- GERAÇÃO DO ARQUIVO ZIP PARA DOWNLOAD ÚNICO ---
                st.success(f"Lote {batch_name} processado com sucesso!")
                
                # Criando um buffer em memória para o arquivo ZIP
                zip_buffer = BytesIO()
                
                # Construindo o ZIP em memória (ZIP_DEFLATED aplica compressão)
                with zipfile.ZipFile(zip_buffer, "w", zipfile.ZIP_DEFLATED) as zip_file:
                    # Escreve o CSV do SISHLA
                    zip_file.writestr(
                        f"resultados_{batch_name}_sishla_format.csv", 
                        df_sishla.to_csv(sep=';', index=False, encoding='utf-8')
                    )
                    
                    # Escreve o CSV do REDOME
                    zip_file.writestr(
                        f"resultados_{batch_name}_redome_format.csv", 
                        df_redome.to_csv(sep=';', index=False, encoding='utf-8')
                    )
                    
                    # Escreve o CSV de Qualidade
                    zip_file.writestr(
                        f"resultados_{batch_name}_informacoes_adicionais.csv", 
                        df_quality_metrics.to_csv(sep=';', index=False, encoding='utf-8')
                    )
                
                # Retorna o ponteiro do buffer para o início para que o Streamlit consiga ler
                zip_buffer.seek(0)

                # Apenas 1 botão centralizado
                st.download_button(
                    label=f"📦 Baixar Pacote de Resultados ({batch_name})",
                    data=zip_buffer,
                    file_name=f"relatorios_completos_{batch_name}.zip",
                    mime="application/zip",
                    use_container_width=True # Deixa o botão mais largo e visível
                )
                    
                st.markdown("---") # Linha separadora para o próximo arquivo
                
            except Exception as e:
                st.error(f"Erro ao processar {batch_name}: {e}")
                
            progress_bar.progress((i + 1) / len(uploaded_files))
            
        status_text.text("Processamento concluído!")

        if 'df_quality_metrics' not in st.session_state:
            st.session_state.df_quality_metrics = None
        
        st.session_state.df_quality_metrics = df_quality_metrics

with tab_dash:
    if st.session_state.df_quality_metrics is not None:
        df_qm = st.session_state.df_quality_metrics
        st.header("Análise de Qualidade da Corrida")
        
        # --- KPIs em Destaque ---
        col_kpi1, col_kpi2, col_kpi3 = st.columns(3)
        total_loci = len(df_qm)
        aprovados = len(df_qm[df_qm['locus_review_status'] == 'Approved'])
        taxa_aprovacao = (aprovados / total_loci) * 100
        
        col_kpi1.metric("Total de Loci Processados", total_loci)
        col_kpi2.metric("Loci Aprovados", aprovados)
        col_kpi3.metric("Taxa de Sucesso", f"{taxa_aprovacao:.1f}%")

        # --- GRÁFICO 1: Profundidade de Leitura (Read Depth) ---
        st.subheader("Distribuição de Profundidade (Exon Min)")
        fig_depth = px.box(df_qm, x="locus_name", y="exon_qm_read_depth_min", 
                           color="locus_name", points="all",
                           title="Profundidade Mínima por Locus (Threshold recomendado > 100)")
        fig_depth.add_hline(y=100, line_dash="dash", line_color="red", annotation_text="Ponto de Corte")
        st.plotly_chart(fig_depth, use_container_width=True)

        # --- GRÁFICO 2: Ruído vs Profundidade ---
        st.subheader("Relação Ruído (Noise) x Profundidade")
        fig_noise = px.scatter(df_qm, x="exon_qm_read_depth_min", y="exon_qm_noise_max",
                               color="locus_review_status", hover_data=['sample_name'],
                               title="Loci Aprovados vs Reprovados pelo Ruído")
        st.plotly_chart(fig_noise, use_container_width=True)

        # --- TABELA DE ALERTAS ---
        st.subheader("⚠️ Alertas: Amostras com Baixa Profundidade (< 50x)")
        alertas = df_qm[df_qm['exon_qm_read_depth_min'].astype(float) < 50]
        if not alertas.empty:
            st.warning(f"Encontrados {len(alertas)} potenciais falhas de cobertura.")
            st.dataframe(alertas[['sample_name', 'locus_name', 'exon_qm_read_depth_min', 'locus_review_status']])
        else:
            st.success("Nenhuma amostra abaixo do threshold crítico de 50x.")
            
    else:
        st.info("Aguardando o processamento de arquivos para gerar o dashboard.")