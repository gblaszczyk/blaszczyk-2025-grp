# Copyright (C) 2025 namemcguffin
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

library(Seurat)
library(slingshot)
library(tidyverse)
library(data.table)
library(patchwork)

rsp <- function(so, cluster.col, start.clus) {
  spt <- slingshot(Embeddings(so, "pca"), so[[]][[cluster.col]], start.clus = start.clus)

  list(
    "curves" = Embeddings(so, "umap") |>
      embedCurves(spt, newDimRed = _) |>
      slingCurves(as.df = TRUE) |>
      as.data.table() |>
      setnames(old = c("Lineage", "Order"), new = c("lineage", "order")) |>
      _[order(order)],
    "pseudotime" = slingPseudotime(spt) |>
      as.data.table(keep.rownames = TRUE) |>
      melt(measure.vars = patterns("Lineage[0-9]+")) |>
      setnames(
        old = c("rn", "variable", "value"),
        new = c("cell", "lineage", "pseudotime")
      ) |>
      _[
        ,
        `:=`(lineage = as.integer(str_extract(lineage, "Lineage([0-9]+)", 1)))
      ],
    "inputs" = list(
      Embeddings(so, "umap"),
      Embeddings(so, "pca"),
      so[[cluster.col]]
    ) |>
      reduce(cbind) |>
      as.data.table(keep.rownames = TRUE) |>
      setnames(old = c("rn"), new = c("cell")) |>
      _[],
    "sds" = spt
  )
}

ptp <- function(
    pseudo.out,
    keep.lineages = NULL) {
  if (!is.null(keep.lineages)) {
    pseudo.out[["pseudotime"]] <- pseudo.out[["pseudotime"]][lineage %in% keep.lineages]
    pseudo.out[["curves"]] <- pseudo.out[["curves"]][lineage %in% keep.lineages]
  }

  pseudo.out[["pseudotime"]][
    pseudo.out[["inputs"]][, .(UMAP_1, UMAP_2, cell)],
    on = .(cell)
  ] |>
    ggplot() +
    aes(
      x = UMAP_1,
      y = UMAP_2
    ) +
    geom_point(
      aes(fill = pseudotime),
      stroke = NA,
      shape = 21
    ) +
    geom_path(
      aes(colour = order),
      data = pseudo.out[["curves"]],
      linewidth = 1
    ) +
    facet_wrap(
      vars(lineage),
      labeller = label_both,
      axes = "all",
      axis.labels = "all"
    ) +
    scale_fill_viridis_c(option = "turbo") +
    scale_colour_viridis_c(option = "rocket") +
    labs(fill = "pseudotime", colour = "order") +
    theme(
      panel.grid = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line()
    )
}
