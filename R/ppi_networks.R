#'Function for accessing the STRINGdb API
#'Adheres to official reccomendations for API use
#' @export
string_api <- function(method, params) {
  
  # The API version is fixed for reproducibility - change if necessary
  
  STRING_URL <- "https://version-12-5.string-db.org"  
  
  req  <- httr2::request(paste0(STRING_URL, "/api/tsv/", method))
  req  <- do.call(httr2::req_body_form, c(list(req), lapply(params, as.character)))
  req  <- httr2::req_error(req, is_error = function(resp) FALSE)
  resp <- httr2::req_perform(req)
  Sys.sleep(1)                                
 
  if (httr2::resp_status(resp) == 404) return(NULL)
  if (httr2::resp_status(resp) >= 400)
    stop("STRING API error: ", httr2::resp_body_string(resp), call. = FALSE)
 
  utils::read.delim(text = httr2::resp_body_string(resp), quote = "",
                    comment.char = "", stringsAsFactors = FALSE)
}

#' Map UniProt accessions to STRING identifiers
#' Returns one row per unique input: uniprot_id, string_id, preferred_name
#' (NA for accessions that could not be mapped).
map_uniprot_to_stringdb <- function(uniprot_ids, species = 9606) {
  
  ids   <- unique(trimws(as.character(unlist(uniprot_ids))))
  ids   <- ids[!is.na(ids) & ids != ""]
  query <- sub("-[0-9]+$", "", ids)                    # drop isoform suffix: P04637-2 -> P04637
 
  hits <- string_api("get_string_ids", list(identifiers = paste(query, collapse = "\r"),
                                            species     = species,
                                            echo_query  = 1,
                                            caller_identity = "R_string_network"))
 
  m <- if (is.null(hits)) rep(NA_integer_, length(ids)) else
       match(seq_along(query) - 1, hits$queryIndex)    # queryIndex is 0-based
  mapping <- data.frame(uniprot_id     = ids,
                        string_id      = if (is.null(hits)) NA_character_ else hits$stringId[m],
                        preferred_name = if (is.null(hits)) NA_character_ else hits$preferredName[m],
                        stringsAsFactors = FALSE)
 
  message(sprintf("Mapped %d of %d accessions.", sum(!is.na(mapping$string_id)), length(ids)))
  
  unmapped <- mapping %>%
    filter(is.na(string_id))%>%
    pull(uniprot_id)
  
  for(protein in unmapped){
    message(glue::glue("{protein} was not mapped to string id."))
  }
  
  mapping
  
}

#' Map gene symbols to STRING identifiers
#' Returns one row per unique input: gene_symbol, string_id, preferred_name
#' (NA for symbols that could not be mapped).
#' @export
map_symbols_to_stringdb <- function(gene_symbols, species = 9606) {
  
  ids <- unique(trimws(as.character(unlist(gene_symbols))))
  ids <- ids[!is.na(ids) & ids != ""]
  
  hits <- string_api("get_string_ids", list(identifiers = paste(ids, collapse = "\r"),
                                            species     = species,
                                            limit       = 1,   # best match per symbol
                                            echo_query  = 1,
                                            caller_identity = "R_string_network"))
  
  m <- if (is.null(hits)) rep(NA_integer_, length(ids)) else
       match(seq_along(ids) - 1, hits$queryIndex)      # queryIndex is 0-based
  mapping <- data.frame(gene_symbol    = ids,
                        string_id      = if (is.null(hits)) NA_character_ else hits$stringId[m],
                        preferred_name = if (is.null(hits)) NA_character_ else hits$preferredName[m],
                        stringsAsFactors = FALSE)
  
  message(sprintf("Mapped %d of %d symbols.", sum(!is.na(mapping$string_id)), length(ids)))
  
  unmapped <- mapping %>%
    filter(is.na(string_id)) %>%
    pull(gene_symbol)
  
  for(symbol in unmapped){
    message(glue::glue("{symbol} was not mapped to string id."))
  }
  
  # symbols resolved through an alias (e.g. an outdated symbol)
  renamed <- mapping %>%
    filter(!is.na(preferred_name), toupper(preferred_name) != toupper(gene_symbol))
  
  for(i in seq_len(nrow(renamed))){
    message(glue::glue("{renamed$gene_symbol[i]} was mapped via alias to {renamed$preferred_name[i]}."))
  }
  
  # several symbols collapsing onto the same STRING protein
  dups <- mapping %>%
    filter(!is.na(string_id), duplicated(string_id) | duplicated(string_id, fromLast = TRUE))
  
  if (nrow(dups) > 0)
    message(glue::glue("These symbols map to the same STRING protein: ",
                       "{paste(dups$gene_symbol, collapse = ', ')}"))
  
  mapping
}

#' Function for building a string network
#' @export
build_string_network <- function(ids, id_type = c("up", "symbol"), species = 9606,
                                 network_type = "physical",
                                 min_score = 0.9, file = NULL, keep.singletons = TRUE) {
  
  id_type <- match.arg(id_type)
  mapping <- switch(id_type,
                    up     = map_uniprot_to_stringdb(ids, species),
                    symbol = map_symbols_to_stringdb(ids, species))
  
  string_ids <- unique(stats::na.omit(mapping$string_id))
  if (length(string_ids) < 2) stop("Need at least 2 mapped proteins.", call. = FALSE)
 
  net <- string_api("network", list(identifiers    = paste(string_ids, collapse = "\r"),
                                    species        = species,
                                    network_type   = network_type,
                                    required_score = round(min_score * 1000),
                                    caller_identity = "R_string_network"))
 
  edges <- if (is.null(net) || nrow(net) == 0) {
    data.frame(node1 = character(0), node2 = character(0), score = numeric(0),
               string_id1 = character(0), string_id2 = character(0))
  } else {
    data.frame(node1      = net$preferredName_A,
               node2      = net$preferredName_B,
               score      = as.numeric(net$score),
               string_id1 = net$stringId_A,
               string_id2 = net$stringId_B,
               stringsAsFactors = FALSE)
  }
 
  # keep each edge once (A-B, not also B-A) and apply the score cut-off
  pair  <- ifelse(edges$string_id1 < edges$string_id2,
                  paste(edges$string_id1, edges$string_id2),
                  paste(edges$string_id2, edges$string_id1))
  edges <- edges[!duplicated(pair) & edges$score >= min_score, ]
 
  # keep proteins without any edge as singletons: a self-edge with score NA
  connected  <- unique(c(edges$string_id1, edges$string_id2))
  singletons <- mapping[!is.na(mapping$string_id) & !mapping$string_id %in% connected, ]
  singletons <- singletons[!duplicated(singletons$string_id), ]
  
  if (nrow(singletons) > 0 && keep.singletons) {
    edges <- rbind(edges, data.frame(node1      = singletons$preferred_name,
                                     node2      = singletons$preferred_name,
                                     score      = NA_real_,
                                     string_id1 = singletons$string_id,
                                     string_id2 = singletons$string_id,
                                     stringsAsFactors = FALSE))
  }
  rownames(edges) <- NULL
 
  message(sprintf("%s network: %d edges with score >= %.2f, plus %d singletons.",
                  network_type, sum(!is.na(edges$score)), min_score, nrow(singletons)))
 
  if (!is.null(file))
    utils::write.table(edges, file, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
 
  attr(edges, "mapping") <- mapping
  edges
}

#'Make a simple network graph with igraph & ggraph
#' @export
plot_string_net_simple <- function(table, c="orange", layout="fr"){
  
  g <- igraph::graph_from_data_frame(
  table,
  directed = FALSE,
  )
    
  ggraph::ggraph(g, layout = layout) +
  ggraph::geom_edge_link(aes(width = score), alpha = 0.3) +
  ggraph::geom_node_point(size = 5, fill=c, color="black", shape=21) +
  ggraph::geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  ggraph::scale_edge_width(range = c(0.2, 2),name = "STRINGdb score") +
  ggraph::theme_graph(base_family = "sans")+
  theme(
    legend.position = "top",
    legend.box      = "vertical",
    legend.box.just = "center"   # optional: align the rows on the left
  )
}


#' Make a simple network graph with igraph & ggraph
#' Nodes are highlighted by "color_by" and colors are set by "cmap"
#' @export
plot_string_net_highlight <- function(table, meta, color_by, name, cmap, layout="fr"){
  g <- igraph::graph_from_data_frame(
  table,
  directed = FALSE,
  vertices = meta
  )
    
  p <- ggraph::ggraph(g, layout = layout) +
  ggraph::geom_edge_link(aes(width = score), alpha = 0.3) +
  ggraph::geom_node_point(aes(fill = !!sym(color_by)), size = 5, color="black", shape=21) +
  ggraph::geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  ggraph::scale_edge_width(range = c(0.2, 2),name = "STRINGdb score") +
  ggraph::theme_graph(base_family = "sans")+
  scale_fill_manual(values=cmap, name=name)+
  theme(
    legend.position = "top",
    legend.box      = "vertical",
    legend.box.just = "center"   # optional: align the rows on the left
  )
  
  return(p)
  
}