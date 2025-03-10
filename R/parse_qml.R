#' Parse QML stylesheet to table of hex values
#' @param qml_path Character. Path of qml stylesheet
#' @details Parse stylesheet
#' @return A data.table containing columns classification and colour
#' @importFrom xml2 read_xml xml_find_all xml_attr
#' @import data.table
#' @export

parse_qml <- function(qml_path){
  doc <- read_xml(qml_path)
  categories <- xml_find_all(doc, './renderer-v2/categories/category')
  labels <- xml_attr(categories,"label")
  symbols <- xml_find_all(doc, './renderer-v2/symbols/symbol/layer/Option/Option')
  temp <- symbols[xml_attr(symbols,"name") == "color"]
  cols <- xml_attr(temp, "value")
  cols <- gsub(",rgb:.*","",cols)
  
  dat <- data.table(classification = labels, fill = cols)
  dat[,c("R","G","B","A") := tstrsplit(fill,",")]
  dat[,colour := rgb(R,G,B,A, maxColorValue = 255)]
  dat <- dat[-1,.(classification,colour)]
  return(dat)
}
