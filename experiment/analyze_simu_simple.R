library(tidyr)
library(dplyr)
library(readr)
library(ggplot2)

n = 10
pq = c(10, 50, 100)
rt = 3
cm = 0

dir_path <- file.path("experiment", "results")

res <- do.call(
  rbind, lapply(
    pq, \(p) {
      file_id <- paste0(
        "_n",n,"_p",p,"_q",p,
        "_rt",rt,
        "_cm",as.integer(cm)
      )
      metafile_name <- paste0("meta-simu-simple", file_id, ".txt")
      file_name <- paste0("res-simu-simple", file_id, ".txt")
      header <- read.table(file.path(dir_path, metafile_name))[['V1']]
      read.csv(file.path(dir_path, file_name), col.names = header, header=FALSE)
    }
  )
)

ggplot(res) +
  geom_point(aes(x = p, y = mseU, color=p))
ggplot(res) +
  geom_point(aes(x = p, y = mseV, color=p))
ggplot(res) +
  geom_point(aes(x = p, y = mseK, color=p))

ggplot(res) +
  geom_point(aes(x = p, y = relmseX, color=p))


