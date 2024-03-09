#### loading required libraries (there might be more libraries loaded than required)
load.lib = c("this.path", "spatstat", "ggplot2")

install.lib = load.lib[!load.lib %in% installed.packages()]
for(lib in install.lib) install.packages(lib, dependencies=TRUE)
sapply(load.lib, require, character=TRUE)

