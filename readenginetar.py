#Input data
import pandas as pd
import glob
from xml_read_engine_functions import *
from mac_request_api import *
from format_dfs_to_import import *

def read_tar_files(input_folder):
  extension = '/*.xml'
  filenames = []
  for filename in glob.glob(input_folder + extension):
    filenames.append(filename)
  return filenames

def extract_metrics_xml(filename):  
  roots = []
  #for filename in filenames:
  tree = ET.parse(filename)  # Read file
  root = tree.getroot()  # Parse XML
  roots.append(root)


    # Concatenating all lines in a DataFrame
  df_cols_id                = ['project', 'software_name','software_version','folder_name', 'fastq_filename', 'sample_name', 'locus_name', 'imgt_version' ]

  df_cols_typing            = ['typing_conclusive', 'typing_result','typing_twofields','typing_threefields','typing_pgroup','nmdp_typing_allele1','nmdp_typing_allele2']

  df_cols_sample_reads      = ['sample_map_accepted', 'sample_map_used', 'sample_map_total','sample_map_percentage']

  df_cols_locus_reads       = ['locus_map_accepted', 'locus_map_used','locus_map_total','locus_map_percentage', 'locus_read_length_average', 'locus_read_depth_analysis_min', 'locus_read_depth_analysis_median']

  df_cols_allele_qm_overall = ['allele_one_ext', 'allele_two_ext', 'allele_one_mm', 'allele_two_mm']
  
  df_review_cols            = ['locus_review_status', 'locus_review_level',
                               'locus_first_reviewer_name', 'locus_first_reviewer_datetime', 'locus_first_reviewer_action',
                               'locus_second_reviewer_name', 'locus_second_reviewer_datetime', 'locus_second_reviewer_action'
                               ]

  df_cols_allele_qm_amp     = ['amp_qm_read_depth_min', 'amp_qm_read_depth_max', 'amp_qm_read_depth_average', 'amp_qm_read_depth_median', 'amp_qm_read_depth_sd',
                              'amp_qm_noise_delta_ston','amp_qm_noise_min', 'amp_qm_noise_max', 'amp_qm_noise_average', 'amp_qm_noise_median', 'amp_qm_noise_std',
                              'amp_qm_secondbase_min', 'amp_qm_secondbase_max', 'amp_qm_secondbase_average', 'amp_qm_secondbase_median', 'amp_qm_secondbase_sd', 'amp_qm_secondbase_count']

  df_cols_allele_qm_exon    = ['exon_qm_read_depth_min', 'exon_qm_read_depth_max', 'exon_qm_read_depth_average', 'exon_qm_read_depth_median', 'exon_qm_read_depth_sd',
                              'exon_qm_noise_delta_ston','exon_qm_noise_min', 'exon_qm_noise_max', 'exon_qm_noise_average', 'exon_qm_noise_median', 'exon_qm_noise_std',
                              'exon_qm_secondbase_min', 'exon_qm_secondbase_max', 'exon_qm_secondbase_average', 'exon_qm_secondbase_median', 'exon_qm_secondbase_sd', 'exon_qm_secondbase_count']

  df_cols_allele_qm_core    = ['core_qm_read_depth_min', 'core_qm_read_depth_max', 'core_qm_read_depth_average', 'core_qm_read_depth_median', 'core_qm_read_depth_sd',
                              'core_qm_noise_delta_ston','core_qm_noise_min', 'core_qm_noise_max', 'core_qm_noise_average', 'core_qm_noise_median', 'core_qm_noise_std',
                              'core_qm_secondbase_min', 'core_qm_secondbase_max', 'core_qm_secondbase_average', 'core_qm_secondbase_median', 'core_qm_secondbase_sd', 'core_qm_secondbase_count']


  df_cols                   = df_cols_id + df_review_cols +df_cols_typing + df_cols_sample_reads + df_cols_locus_reads + df_cols_allele_qm_amp + df_cols_allele_qm_exon + df_cols_allele_qm_core + df_cols_allele_qm_overall


  rows = []
  for root in roots:
    project = xml_text_value(root, 'ProjectName' )
    software_name = xml_text_value(root, 'AnalysisSoftware/Software' )
    software_version = xml_text_value(root, 'AnalysisSoftware/Version' )

    samples = xml_many_elements(root, './Samples/Sample')

    for sample in samples:
      xml                 = sample
      path                = 'Name'
      fastq_filename      = xml_text_value(xml, path)
      sample_name         = fastq_filename.split('-')[0]
      path                = "Mappability"
      map_attrib          = ['AcceptedReads', 'UsedReads', 'TotalReads', 'PercentageUsed' ]
      sample_map_metrics  = many_xml_attrib( xml, path, map_attrib)
      sample_map_accepted, sample_map_used, sample_map_total, sample_map_percentage = sample_map_metrics

      path                = "Files/File"
      folder_name         = xml_many_text_elements(xml, path)[0].split('\\')[-2]

      path              = './Loci/Locus'
      loci = xml_many_elements(xml, path)

      for locus in loci:
        xml               = locus

        path              = 'Name'
        locus_name        = xml_text_value(xml, path )

        #HLA Locus IMGT version database
        path              = 'AlleleDB/Version'
        imgt_version      = xml_text_value(xml, path )        

        #HLA Locus quality metrics
        path              = "Mappability"
        locus_map_metrics = many_xml_attrib( xml, path, map_attrib)
        locus_map_accepted, locus_map_used, locus_map_total, locus_map_percentage = locus_map_metrics

        path              = "ReadLength"
        locus_read_length_average = xml_one_element(xml, path, "Mean" )

        statistic_metrics = ["Minimum", "Median"]
        path              = "ReadDepth"
        locus_sta_metrics = many_xml_attrib( xml, path, statistic_metrics)
        locus_read_depth_analysis_min, locus_read_depth_analysis_median = locus_sta_metrics

        # #HLA Locus Module Review - start
        path              = "ReviewList"
        locus_review_status = xml_one_element(xml, path, "CurrentApprovalStatus" )
        locus_review_level  = xml_one_element(xml, path, "Level" )

        path              = "ReviewList/Review"
        reviews           = xml_many_elements(xml, path)

        locus_first_reviewer_name, locus_first_reviewer_datetime, locus_first_reviewer_action = "", "", ""
        locus_second_reviewer_name, locus_second_reviewer_datetime, locus_second_reviewer_action = "", "", ""

        for review in reviews:
          
          reviewer_action = review.attrib.get("Action")

          path              = "User"
          reviewer_name     = xml_one_element(review, path, "Name" )
          reviewer_datetime = xml_one_element(review, path, "DateTime" )
          reviewer_level     = xml_one_element(review, path, "Level" )

          if reviewer_level == "First reviewer":
            locus_first_reviewer_name     = reviewer_name
            locus_first_reviewer_datetime = reviewer_datetime
            locus_first_reviewer_action   = reviewer_action
          else:
            locus_second_reviewer_name     = reviewer_name
            locus_second_reviewer_datetime = reviewer_datetime
            locus_second_reviewer_action   = reviewer_action
        
        # #HLA Locus Module Review - END

        typing_paths      = ["ConclusiveTypingResult/Typing/GLString",
                            "TypingResult/GLString",
                            "TypingResultTwoFields/GLString",
                            "TypingResultThreeFields/GLString",
                            "TypingResultPGroup/GLString"]
        locus_typings     = []

        for path in typing_paths:
          typing = xml_text_value(xml, path)
          locus_typings.append(typing)

        typing_conclusive, typing_result, typing_twofields, typing_threefields, typing_pgroup = locus_typings

        nmdp_paths    = ["TypingResult/NMDPCodes/allele_1",
                        "TypingResult/NMDPCodes/allele_2"]
        attr          = 'Name'
        nmdp_typings  = []

        for path in nmdp_paths:
          typing = xml_one_element(xml, path, attr)
          nmdp_typings.append(typing)

        nmdp_allele1, nmdp_allele2 = nmdp_typings

        # Allele quality metrics overall
        path          = "./AnnotatedAlleles/AnnotatedAllele/Extended"
        allele_one_ext,allele_two_ext = xml_many_text_elements( locus, path)

        path               =  "./AnnotatedAlleles/AnnotatedAllele/Mismatches"
        allele_one_mm, allele_two_mm = xml_many_text_elements(locus, path )
        
        # Inicialize todas as variáveis de métricas de região antes do loop
        amp_qm_noise_delta_ston, amp_qm_read_depth_min, amp_qm_read_depth_max, amp_qm_read_depth_average, amp_qm_read_depth_median, amp_qm_read_depth_sd = [None] * 6
        amp_qm_noise_min, amp_qm_noise_max, amp_qm_noise_average, amp_qm_noise_median, amp_qm_noise_std = [None] * 5
        amp_qm_secondbase_min, amp_qm_secondbase_max, amp_qm_secondbase_average, amp_qm_secondbase_median, amp_qm_secondbase_sd, amp_qm_secondbase_count = [None] * 6

        exon_qm_noise_delta_ston, exon_qm_read_depth_min, exon_qm_read_depth_max, exon_qm_read_depth_average, exon_qm_read_depth_median, exon_qm_read_depth_sd = [None] * 6
        exon_qm_noise_min, exon_qm_noise_max, exon_qm_noise_average, exon_qm_noise_median, exon_qm_noise_std = [None] * 5
        exon_qm_secondbase_min, exon_qm_secondbase_max, exon_qm_secondbase_average, exon_qm_secondbase_median, exon_qm_secondbase_sd, exon_qm_secondbase_count = [None] * 6

        core_qm_noise_delta_ston, core_qm_read_depth_min, core_qm_read_depth_max, core_qm_read_depth_average, core_qm_read_depth_median, core_qm_read_depth_sd = [None] * 6
        core_qm_noise_min, core_qm_noise_max, core_qm_noise_average, core_qm_noise_median, core_qm_noise_std = [None] * 5
        core_qm_secondbase_min, core_qm_secondbase_max, core_qm_secondbase_average, core_qm_secondbase_median, core_qm_secondbase_sd, core_qm_secondbase_count = [None] * 6
        
        path                =   "./PreferredAnalysisRegionsQualityMetrics/RegionQualityMetrics"
        qm_regions          =   xml_many_elements(xml, path)

## Region quality metrics - start###########################
        for qm in qm_regions:
          xml     = qm
          region  = qm.attrib.get('region')

          attr    = 'minimumHeterozygousVsMaxTotalNoise'
          path    = './Noise'
          delta   = xml_one_element(xml, path, attr)

          attrs   = ['minimum', 'maximum', 'average', 'median', 'sd']

          path    = './InputData/ReadDepth'
          region_read_depth_metrics = []

          for attr in attrs:
            metric_data = xml_one_element(xml, path, attr)
            region_read_depth_metrics.append(metric_data)

          path    = './Noise/TotalNoiseLevel'
          noise_metrics = []
          for attr in attrs:
            metric_data = xml_one_element(xml, path, attr)
            noise_metrics.append(metric_data)


          attrs       = ['minimum', 'maximum', 'average', 'median', 'sd', 'count']
          path        = './Noise/SecondBaseStatisticsAtHeterozygotes'
          secondbase_metrics = []
          for attr in attrs:
            metric_data = xml_one_element(xml, path, attr)
            secondbase_metrics.append(metric_data)

          if region == 'Amplicon':
            amp_qm_noise_delta_ston = delta
            amp_qm_read_depth_min, amp_qm_read_depth_max, amp_qm_read_depth_average, amp_qm_read_depth_median, amp_qm_read_depth_sd = region_read_depth_metrics
            amp_qm_noise_min, amp_qm_noise_max, amp_qm_noise_average, amp_qm_noise_median, amp_qm_noise_std = noise_metrics
            amp_qm_secondbase_min, amp_qm_secondbase_max, amp_qm_secondbase_average, amp_qm_secondbase_median, amp_qm_secondbase_sd, amp_qm_secondbase_count = secondbase_metrics

          if region == 'Exon+':
            exon_qm_noise_delta_ston  = delta
            exon_qm_read_depth_min, exon_qm_read_depth_max, exon_qm_read_depth_average, exon_qm_read_depth_median, exon_qm_read_depth_sd = region_read_depth_metrics
            exon_qm_noise_min, exon_qm_noise_max, exon_qm_noise_average, exon_qm_noise_median, exon_qm_noise_std = noise_metrics
            exon_qm_secondbase_min, exon_qm_secondbase_max, exon_qm_secondbase_average, exon_qm_secondbase_median, exon_qm_secondbase_sd, exon_qm_secondbase_count = secondbase_metrics

          if region == 'Core+':
            core_qm_noise_delta_ston  = delta
            core_qm_read_depth_min, core_qm_read_depth_max, core_qm_read_depth_average, core_qm_read_depth_median, core_qm_read_depth_sd = region_read_depth_metrics
            core_qm_noise_min, core_qm_noise_max, core_qm_noise_average, core_qm_noise_median, core_qm_noise_std = noise_metrics
            core_qm_secondbase_min, core_qm_secondbase_max, core_qm_secondbase_average, core_qm_secondbase_median, core_qm_secondbase_sd, core_qm_secondbase_count = secondbase_metrics

        rows.append({
                  'project': project,
                  'folder_name':folder_name,
                  'fastq_filename': fastq_filename,
                  'sample_name': sample_name,
                  'software_name': software_name,
                  'software_version': software_version,
                  'imgt_version': imgt_version,

                  # review
                  "locus_review_status":locus_review_status,
                  "locus_review_level":locus_review_level,
                  "locus_first_reviewer_name":locus_first_reviewer_name,
                  "locus_first_reviewer_datetime":locus_first_reviewer_datetime,
                  "locus_first_reviewer_action":locus_first_reviewer_action,

                  "locus_second_reviewer_name":locus_second_reviewer_name,
                  "locus_second_reviewer_datetime":locus_second_reviewer_datetime,
                  "locus_second_reviewer_action":locus_second_reviewer_action,

                  # Locus Typing
                  'typing_conclusive': typing_conclusive,
                  'typing_result': typing_result,
                  'typing_twofields': typing_twofields,
                  'typing_threefields': typing_threefields,
                  'typing_pgroup': typing_pgroup,
                  'nmdp_typing_allele1': nmdp_allele1,
                  'nmdp_typing_allele2': nmdp_allele2,
                  # Sample Metrics
                  'sample_map_accepted': sample_map_accepted,
                  'sample_map_used': sample_map_used,
                  'sample_map_total':sample_map_total,
                  'sample_map_percentage':sample_map_percentage,

                  #HLA Locus
                  'locus_name': locus_name,
                  # Locus Metrics
                  'locus_map_accepted': locus_map_accepted,
                  'locus_map_used': locus_map_used,
                  'locus_map_total': locus_map_total,
                  'locus_map_percentage': locus_map_percentage,
                  'locus_read_length_average': locus_read_length_average,
                  'locus_read_depth_analysis_min': locus_read_depth_analysis_min,
                  'locus_read_depth_analysis_median': locus_read_depth_analysis_median,
                  'review_status': locus_review_status,

                  # Allele quality metrics overall
                  'allele_one_ext': allele_one_ext,
                  'allele_two_ext': allele_two_ext,
                  'allele_one_mm': allele_one_mm,
                  'allele_two_mm': allele_two_mm,

                  #  Allele quality metrics by region
                  # AMPLICON
                  'amp_qm_read_depth_min': amp_qm_read_depth_min,
                  'amp_qm_read_depth_max': amp_qm_read_depth_max,
                  'amp_qm_read_depth_average': amp_qm_read_depth_average,
                  'amp_qm_read_depth_median': amp_qm_read_depth_median,
                  'amp_qm_read_depth_sd': amp_qm_read_depth_sd,

                  'amp_qm_noise_delta_ston': amp_qm_noise_delta_ston,
                  'amp_qm_noise_min': amp_qm_noise_min,
                  'amp_qm_noise_max': amp_qm_noise_max,
                  'amp_qm_noise_average': amp_qm_noise_average,
                  'amp_qm_noise_median': amp_qm_noise_median,
                  'amp_qm_noise_std': amp_qm_noise_std,
                  'amp_qm_secondbase_min': amp_qm_secondbase_min,
                  'amp_qm_secondbase_max': amp_qm_secondbase_max,
                  'amp_qm_secondbase_average': amp_qm_secondbase_average,
                  'amp_qm_secondbase_median': amp_qm_secondbase_median,
                  'amp_qm_secondbase_sd': amp_qm_secondbase_sd,
                  'amp_qm_secondbase_count': amp_qm_secondbase_count,
                  # EXON
                  'exon_qm_read_depth_min': exon_qm_read_depth_min,
                  'exon_qm_read_depth_max': exon_qm_read_depth_max,
                  'exon_qm_read_depth_average': exon_qm_read_depth_average,
                  'exon_qm_read_depth_median': exon_qm_read_depth_median,
                  'exon_qm_read_depth_sd': exon_qm_read_depth_sd,

                  'exon_qm_noise_delta_ston': exon_qm_noise_delta_ston,
                  'exon_qm_noise_min': exon_qm_noise_min,
                  'exon_qm_noise_max': exon_qm_noise_max,
                  'exon_qm_noise_average': exon_qm_noise_average,
                  'exon_qm_noise_median': exon_qm_noise_median,
                  'exon_qm_noise_std': exon_qm_noise_std,
                  'exon_qm_secondbase_min': exon_qm_secondbase_min,
                  'exon_qm_secondbase_max': exon_qm_secondbase_max,
                  'exon_qm_secondbase_average': exon_qm_secondbase_average,
                  'exon_qm_secondbase_median': exon_qm_secondbase_median,
                  'exon_qm_secondbase_sd': exon_qm_secondbase_sd,
                  'exon_qm_secondbase_count': exon_qm_secondbase_count,

                  # CORE
                  'core_qm_read_depth_min': core_qm_read_depth_min,
                  'core_qm_read_depth_max': core_qm_read_depth_max,
                  'core_qm_read_depth_average': core_qm_read_depth_average,
                  'core_qm_read_depth_median': core_qm_read_depth_median,
                  'core_qm_read_depth_sd': core_qm_read_depth_sd,

                  'core_qm_noise_delta_ston': core_qm_noise_delta_ston,
                  'core_qm_noise_min': core_qm_noise_min,
                  'core_qm_noise_max': core_qm_noise_max,
                  'core_qm_noise_average': core_qm_noise_average,
                  'core_qm_noise_median': core_qm_noise_median,
                  'core_qm_noise_std': core_qm_noise_std,
                  'core_qm_secondbase_min': core_qm_secondbase_min,
                  'core_qm_secondbase_max': core_qm_secondbase_max,
                  'core_qm_secondbase_average': core_qm_secondbase_average,
                  'core_qm_secondbase_median': core_qm_secondbase_median,
                  'core_qm_secondbase_sd': core_qm_secondbase_sd,
                  'core_qm_secondbase_count': core_qm_secondbase_count
                  })


  df = pd.DataFrame(rows, columns = df_cols)
  return df

