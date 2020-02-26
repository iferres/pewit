#' @name PewitR6
#' @importFrom R6 R6Class
#' @import pagoo
#' @importFrom gggenes geom_gene_arrow theme_genes
#' @importFrom S4Vectors List
#' @importFrom ggplot2 ggplot aes theme element_blank
#' @importFrom gggenes geom_gene_arrow theme_genes
PewitR6 <- R6Class('PewitR6',

                   inherit = pagoo::PgR6MS,

                   public = list(

                     initialize = function(data,
                                           org_meta,
                                           cluster_meta,
                                           core_level = 95,
                                           sep = '__',
                                           verbose = TRUE){

                       super$initialize(data = data,
                                        org_meta,
                                        cluster_meta,
                                        core_level = core_level,
                                        sep = sep,
                                        verbose = verbose)

                     },

                     gene_context = function(cluster,
                                             upstream = 2,
                                             downstram = 2,
                                             orient = TRUE){

                       # Compute
                       group <- self$gene[[cluster]]
                       group$org <- sub("__.+", "", group$gid)
                       data <- private$.data
                       ap <- apply(group, 1, function(x) {
                         contig <- data[which(data$org == x[["org"]] &
                                                data$contig == x[["contig"]]), ,drop = FALSE]
                         strand <- x[["strand"]]
                         gix <- which(contig$gid == x[["gid"]])

                         if (strand == "+"){
                           gixup <- gix - upstream
                           gixup <- ifelse(gixup < 1, 1, gixup)
                           gixdown <- gix + downstram
                           gixdown <- ifelse (gixdown > dim(contig)[1], dim(contig)[1], gixdown)
                           contig2 <- contig[gixup:gixdown, ,drop=F]
                         }else{
                           gixup <- gix - downstream
                           gixup <- ifelse(gixup < 1, 1, gixup)
                           gixdown <- gix + upstream
                           gixdown <- ifelse (gixdown > dim(contig)[1], dim(contig)[1], gixdown)
                           contig2 <- contig[gixup:gixdown, ,drop=F]
                         }

                         # orient
                         if (orient & strand == "-"){
                           contig2$from <- contig2$from * -1
                           contig2$to <- contig2$to * -1
                         }

                         # center
                         if (strand == "+"){
                           dif <- as.integer(contig2$from[which(contig2$gid == x[["gid"]])])
                         }else{
                           dif <- as.integer(contig2$to[which(contig2$gid == x[["gid"]])])
                         }
                         contig2$from <- contig2$from - dif
                         contig2$to <- contig2$to - dif

                         return(contig2)
                       })

                       df <- unlist(List(ap))
                       df

                     },

                     gg_gene_context = function(cluster,
                                                upstram = 2,
                                                downstram = 2,
                                                fill = "cluster",
                                                orient = TRUE){

                       df <- as.data.frame(self$gene_context(cluster, upstram, downstram, orient))
                       df$direction <- ifelse(df$strand == "+", 1, -1)
                       ggplot2(df, aes(xmin = from, xmax = to, y = org, fill = cluster, forward = direction)) +
                         geom_gene_arrow() +
                         theme_genes() +
                         theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

                     }

                   )


)
