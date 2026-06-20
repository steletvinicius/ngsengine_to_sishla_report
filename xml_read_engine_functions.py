# Read XML data
import xml.etree.ElementTree as ET

def xml_element_valid (xml_element, path ):
  return True if xml_element.find(path) is not None else False

def many_xml_attrib( xml_element, path, attributes):
  results = []
  for attr in attributes:
    data = xml_one_element( xml_element, path, attr )
    results.append(data)
  return results

def xml_one_element( xml_element, path, attrib):
  if xml_element_valid( xml_element, path):
    return xml_element.find(path).attrib.get(attrib)
  else:
    return None

def xml_many_elements( xml_element, path ):
  return xml_element.findall(path)

def xml_text_value( xml_element, path  ):
  if xml_element_valid(xml_element, path):
    return xml_element.find(path).text
  else:
    return None

def xml_many_text_elements(xml_element, path):
  list_elements = xml_many_elements( xml_element, path)
  if len(list_elements) > 1:
    text_one  = list_elements[0].text
    text_two  = list_elements[1].text
  elif len(list_elements) == 1:
    text_one  = list_elements[0].text
    text_two  = 'NA'
  else:
    text_one  = 'NA'
    text_two  = 'NA'
  return(text_one, text_two)
