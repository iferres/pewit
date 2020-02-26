#' @name PewitR6
#' @importFrom R6 R6Class
#' @import pagoo
#' @importFrom gggenes geom_gene_arrow theme_genes
#' @importFrom S4Vectors List
#' @importFrom ggplot2 ggplot aes theme guides guide_legend
#' @importFrom gggenes geom_gene_arrow theme_genes
PewitR6 <- R6Class('PewitR6',

                   inherit = pagoo::PgR6MS,

                   public = list(

                     initialize = function(data,
                                           org_meta,
                                           cluster_meta,
                                           core_level = 95,
                                           sep = '__',
                                           sequences,
                                           verbose = TRUE){

                       super$initialize(data = data,
                                        org_meta,
                                        cluster_meta,
                                        core_level = core_level,
                                        sep = sep,
                                        sequences = sequences,
                                        verbose = verbose)

                     },

                     gene_context = function(cluster,
                                             upstream = 2,
                                             downstream = 2){

                       # Compute
                       group <- self$genes[[cluster]]
                       sep <- private$.sep
                       group$org <- sub(paste0(sep, ".+"), "", group$gid)
                       data <- private$.data
                       ap <- apply(group, 1, function(x) {
                         contig <- data[which(data$org == x[["org"]] &
                                                data$contig == x[["contig"]]), ,drop = FALSE]
                         strand <- x[["strand"]]
                         gix <- which(contig$gid == x[["gid"]])

                         if (strand == "+"){
                           gixup <- gix - upstream
                           gixup <- ifelse(gixup < 1, 1, gixup)
                           gixdown <- gix + downstream
                           gixdown <- ifelse (gixdown > dim(contig)[1], dim(contig)[1], gixdown)
                           contig2 <- contig[gixup:gixdown, ,drop=F]
                         }else{
                           gixup <- gix - downstream
                           gixup <- ifelse(gixup < 1, 1, gixup)
                           gixdown <- gix + upstream
                           gixdown <- ifelse (gixdown > dim(contig)[1], dim(contig)[1], gixdown)
                           contig2 <- contig[gixup:gixdown, ,drop=F]
                         }

                         return(contig2)
                       })

                       lst <- List(ap)
                       ulst <- unlist(lst)
                       return(split(ulst, as.character(ulst$org)))

                     },

                     gg_gene_context = function(cluster,
                                                upstream = 2,
                                                downstream = 2,
                                                fill = "cluster",
                                                orient = TRUE){

                       group <- self$genes[[cluster]]
                       strands <- group$strand
                       gids <- group$gid
                       contexts <- self$gene_context(cluster, upstream, downstream)

                       mp <- mapply(function(context, strand, gid, orient) {
                         x2 <- context

                         # orient
                         if (orient & strand == "-"){
                           x2$from <- context$from * -1
                           x2$to <- context$to * -1
                         }

                         # center
                         if (strand == "+"){
                           dif <- as.integer(x2$from[which(x2$gid == gid)])
                         }else{
                           dif <- as.integer(x2$to[which(x2$gid == gid)])
                         }
                         x2$from <- x2$from - dif
                         x2$to <- x2$to - dif

                         return(x2)

                       },
                       context = contexts,
                       strand = strands,
                       gid = gids,
                       MoreArgs = list(orient = orient))

                       df <- as.data.frame(unlist(List(mp)))
                       df$geneName <- sub("^$", NA_character_, df$geneName)
                       df$Pfam_Arch <- sub("^NOARCH_\\d+", NA_character_, df$Pfam_Arch)
                       df$direction <- ifelse(df$strand == "+", 1, -1)

                       ggplot(df, aes(xmin = from, xmax = to, y = org, fill = df[[fill]], forward = direction)) +
                         geom_gene_arrow() +
                         theme_genes() +
                         guides(fill = guide_legend(title = fill))

                     }

                   )


)
