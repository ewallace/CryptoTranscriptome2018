# sharedScripts.R
# common packages and functions used in CryptoTranscriptome2018

library(tidyverse)
library(cowplot)
library(ggrepel)
library(ggseqlogo)
library(Biostrings)

## Graphics settings & functions

theme_set(theme_cowplot(font_size=12) + 
              theme(strip.background = element_blank()) )

geom_diagline <- function(linetype='solid',size=0.1,colour="grey20",...) {
    geom_abline(slope=1,intercept=0,linetype=linetype,colour=colour)
}

scientific_10 <- function(x) {
    xout <- gsub("1e", "10^{", format(x),fixed=TRUE)
    xout <- gsub("{-0", "{-", xout,fixed=TRUE)
    xout <- gsub("{+", "{", xout,fixed=TRUE)
    xout <- gsub("{0", "{", xout,fixed=TRUE)
    xout <- paste(xout,"}",sep="")
    return(parse(text=xout))
}

scale_x_log10nice <- function(name=waiver(),omag=seq(-10,20),...) {
    breaks10 <- 10^omag
    scale_x_log10(name,breaks=breaks10,labels=scientific_10(breaks10),...)
}

scale_y_log10nice <- function(name=waiver(),omag=seq(-10,20),...) {
    breaks10 <- 10^omag
    scale_y_log10(name,breaks=breaks10,labels=scientific_10(breaks10),...)
}

scale_loglog <- function(...) {
    list(scale_x_log10nice(...),scale_y_log10nice(...))
}

## Numeric convenience functions

round3 <- function(x) round(x,digits=3)

logmean <- function(x) {
    lx <- log(x)
    exp( mean( lx[is.finite(lx)] ) )
}

## ATG context functions

read_ATGcontext_table <- function(file) {
    read_tsv(file,comment="#",skip=1,
             col_names=c("Gene","aATG.context","aATG.pos",
                         "d1.context","d1.posTSS","d1.posATG","d1.frame",
                         "d2.context","d2.posTSS","d2.posATG","d2.frame",
                         "u1.context","u1.posTSS","u1.posATG","u1.frame",
                         "u2.context","u2.posTSS","u2.posATG","u2.frame")
    ) %>% 
        filter(!is.na(aATG.context) ) %>%
        mutate(aATG.context=toupper(aATG.context))
}

## Information content functions

infon1 <- function(counts,ssize=4,pseudocount=0.5) {
    # calculate information content for one position
    # as described in https://en.wikipedia.org/wiki/Sequence_logo
    # counts: position frequency matrix
    # ssize: alphabet size
    # pseudocount: pseudocount to add to remove zeroes
    
    counts <- counts + pseudocount # add pseudocount
    ntot = sum(counts)[1] # total number of sequences
    freqs = counts / ntot # freqency
    Hentropy = - sum( freqs * log2(freqs) ) # Shannon entropy 
    ecorrection = ( ssize - 1 ) / ( log(2) *  ntot ) # small-sample correction
    information = log2(ssize) - Hentropy # information content
    return(information)
}

infonall <- function(pfm,ssize=4,pseudocount=0.5) {
    apply(pfm,2,infon1)
}

infontot <- function(pfm,ssize=4,pseudocount=0.5) {
    apply(pfm,2,infon1) %>% sum()
}

infonfromseqs <- function(seqs,sstart=1L,send=-1L,ssize=4,pseudocount=0.5) {
    seqs %>% 
    str_sub(start=sstart,end=send) %>%
    str_to_upper %>%
    consensusMatrix %>%
    infontot(ssize=ssize,pseudocount=pseudocount)
}


PWMscore <- function(seqs,pwm,startl) {
    # calculate scores for a character vector seqs against
    # pwm from given starting position
    # after checking strings are non-missing and start with ACTG
    vapply(str_sub(seqs,start=startl),
           function(ss) {
               if ( is.na(ss) ) {
                   return(NA)
               } else if ( str_detect(ss,"\\A[ACTG]") ) {
                   return( PWMscoreStartingAt(pwm,ss) )
               } else {
                   return(NA) 
               } 
           } ,
           FUN.VALUE = 0)
}

PWMscoren <- function(seqs,pwm=kozak_n_PWM_Sc,startl=8) {
    PWMscore(seqs,pwm,startl) 
}

PWMscorew <- function(seqs,pwm=kozak_w_PWM_Sc,startl=4) {
    PWMscore(seqs,pwm,startl) 
}

is.context <- function(x) str_detect(x,"context")

