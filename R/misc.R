## Software License Agreement (BSD License)
##
## Copyright (c) 2014, Ross Linscott (rossklin@gmail.com)
##
## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are
## met:
##
##     Redistributions of source code must retain the above copyright
##     notice, this list of conditions and the following disclaimer.
##
##     Redistributions in binary form must reproduce the above copyright
##     notice, this list of conditions and the following disclaimer in
##     the documentation and/or other materials provided with the
##     distribution.
##
##     The names of its contributors may not be used to endorse or promote products
##     derived from this software without specific prior written permission.
##
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
## "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
## LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
## A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
## HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
## SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
## LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
## DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
## THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
## (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
## OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

standardise <- function(xs) (xs - mean(na.omit(xs)))/sd(na.omit(xs))
norm2 <- function(xs) norm(as.matrix(xs), type="2")
rescale <- function(xs, scale=1) scale*(xs/norm2(xs))

subsets <- function(n) {
    ss <- expand.grid(rep(list(c(FALSE, TRUE)), n))
    sizes <- rowSums(ss)
    ss[do.call("order", c(list(s=sizes), ss)),]
}

load.as.tt <- function(sel){
  ldf <- load.raw.data(sel)
  tt <- as.time.table(ldf$data, ldf$id.cols, ldf$time.cols, ldf$data.cols)
}

load.raw.data <- function(sel){

indicators <- as.data.frame(matrix(byrow=T, ncol=2,
  c(    "SP.DYN.TFRT.IN", "Births" #Fertility rate (birth per woman)
   ,       "SP.POP.GROW", "Population.Growth" #(annual %)
   ,    "SL.UEM.TOTL.ZS", "Unemployment" #, total (% of total labor force)
   , "SP.POP.65UP.TO.ZS", "Elderly" #Population ages 65 and above (% of total)
   ,    "SP.RUR.TOTL.ZS", "Rural.Population" #(% of total population)
   ,    "AG.LND.FRST.ZS", "Forest.Area" #(% of land area)
   ,    "AG.LND.AGRI.ZS", "Agricultural.Land" #(% of land area)
   ,    "AG.LND.TOTL.K2", "Land.Area" #(sq. km)
   ,       "SP.POP.TOTL", "Population" # (Total)
   )), stringsAsFactors=FALSE)
colnames(indicators) <- c("Code", "Name")

#dataset.raw <- WDI(country="all", indicator=indicators$Code, start=1961, end=2012)
dataset.raw <- read.csv("wdi_cached.csv")
dataset.raw <- setNames(dataset.raw, c("index", "iso2c", "country", "year", indicators$Name))
dataset.raw$iso2c <- as.character(dataset.raw$iso2c)
fix <- is.na(dataset.raw$iso2c)
dataset.raw$iso2c[fix] <- "NA"

if (sel == "all"){
  sel <- indicators$Name
}

res <- NULL
res$id.cols <- "iso2c"
res$time.cols <- "year"
res$data.cols <- sel
res$data <- dataset.raw[,c(res$id.cols, res$time.cols, res$data.cols)]
res$data[,res$id.cols] <- as.character(res$data[,res$id.cols])
res
}