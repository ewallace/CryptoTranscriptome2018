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

## string and codon functions
str_dice <- function(string, width=1L) {
    # Dice a character vector into pieces of fixed length.
    # Useful for breaking a coding sequence of DNA into codons
    starts <- seq(1, nchar(string), width)
    str_sub(string, start = starts, end = starts+width-1)
}


trim_at_last_stop_codon <- function(codonseq, stopcods=c("TGA","TAG","TAA"),re) {
    ncodons <- length(codonseq)
    last_stop <- which(codonseq %in% stopcods) %>% max
    if (!is.finite(last_stop)) {
        return(codonseq) 
    } else if (last_stop == ncodons) {
        return("") 
    } else if (last_stop < ncodons) {
        return( codonseq[(last_stop+1):ncodons] )
    }
}

trim_to_first_nrstart <- function(codonseq, 
                               startcods=c("ATG",
                                           "TTG","CTG","ACG",
                                           "ATT","GTG","ATA")) {
    ncodons <- length(codonseq)
    first_start <- which(codonseq %in% startcods) %>% min
    if (!is.finite(first_start)) {
        return("") 
    } else if (first_start >= 1) {
        return( codonseq[first_start:ncodons] )
    }
}

estimate_3prime_nrstart_CDS <- function(dseq,
                                        stopcods=c("TGA","TAG","TAA"),
                                        startcods=c("ATG",
                                           "TTG","CTG","ACG",
                                           "ATT","GTG","ATA"),
                                        frameright=TRUE) {
    if (frameright) {
        # align frame to right/3' end
        dslength <- nchar(dseq)
        nonframe <- dslength %% 3
        dseq <- str_sub(dseq,start=nonframe+1)
    }
    dseq %>%
        str_dice(width=3L) %>% # dice into codons
        trim_at_last_stop_codon(stopcods=stopcods) %>% # trim 5' end at last stop
        trim_to_first_nrstart(startcods=startcods) %>% # trim 5' end from first start
        paste0(collapse="") # paste into letters
}


threeprime_nrstarts <- function(x,
                                stopcods=c("TGA","TAG","TAA"),
                                startcods=c("ATG",
                                            "TTG","CTG","ACG",
                                            "ATT","GTG","ATA")) {
    as.character(x) %>%
        vapply(estimate_3prime_nrstart_CDS,
               stopcods=stopcods,startcods=startcods,
               FUN.VALUE = "ATG") %>%
        DNAStringSet %>%
        setNames(names(x))
}

# Test these functions
# ORF <- "ATGTTTGGGTAG"
# notORF <- "TGATAATAGAAACCCTTTGGGTAAATG"
# nostop <- "ATGTTTGGGCTC"
# stopstart <- "AAATAGAAAATTGAA"
# startstop <- "ATGTAAAAA"
# startstart <- "ATACCCACG"
# startstopstart <- "ACGAAATGAAAATTGAAA"
# badlength <- "ATGATGAT"
# ORF %>% str_dice(width=3) %>% trim_at_last_stop_codon
# notORF %>% str_dice(width=3) %>% trim_at_last_stop_codon
# nostop %>% str_dice(width=3) %>% trim_at_last_stop_codon
# ORF %>% str_dice(width=3) %>% trim_to_first_nrstart
# notORF %>% str_dice(width=3) %>% trim_to_first_nrstart
# nostop %>% str_dice(width=3) %>% trim_to_first_nrstart
# stopstart %>% str_dice(width=3) %>% trim_to_first_nrstart
# 
# ORF %>% estimate_3prime_nrstart_CDS
# notORF %>% estimate_3prime_nrstart_CDS
# nostop %>% estimate_3prime_nrstart_CDS
# stopstart %>% estimate_3prime_nrstart_CDS
# startstop %>% estimate_3prime_nrstart_CDS
# startstart %>% estimate_3prime_nrstart_CDS
# startstopstart %>% estimate_3prime_nrstart_CDS

expand_threeprime_nrstarts <- function(txinfile,CDSoutfile,protoutfile=NULL,
                                       leftpad=120,rightpad=120,
                                stopcods=c("TGA","TAG","TAA"),
                                startcods=c("ATG",
                                            "TTG","CTG","ACG",
                                            "ATT","GTG","ATA"),
                                trimstop = 1) {
    
    txin <- readDNAStringSet(txinfile)
    fiveUTR <- subseq(txin,start=1,end=leftpad)
    CDS <- subseq(txin,start=(leftpad+1),end=-(rightpad+1))
    
    nrstartfive <- as.character(fiveUTR) %>%
        vapply(estimate_3prime_nrstart_CDS,
               stopcods=stopcods,startcods=startcods,
               FUN.VALUE = "ATG") %>%
        DNAStringSet
    
    expanded_CDS <- xscat(nrstartfive,CDS) %>%
        setNames( names(CDS) %>% 
                      str_extract(boundary("word") ) %>%
                      paste("Ntermext") ) 
    
    writeXStringSet(expanded_CDS,filepath = CDSoutfile)
    
    if (!is.null(protoutfile)) { 
        expanded_CDS %>%
            translate(if.fuzzy.codon = "X") %>%
            subseq(end= -(1+trimstop) ) %>%
            writeXStringSet(filepath = protoutfile)
    }
    
    return(expanded_CDS)
}

## Run mitofates from an R session, after trimming the strings appropriately.

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

## Read mitofates output

read_mitofates <- function(mffile,              
                           mitofates_colnames = 
                               c("Gene", "Prob_preseq", "Pred_preseq", 
                                 "CSite_enzyme", "Net_charge", "Pos_TOM20", "Pos_afhelix", 
                                 "BHHPPP","BPHBHH", "HBHHBb", "HBHHbB", 
                                 "HHBHHB", "HHBPHB", "HHBPHH", "HHBPHP",
                                 "HHHBBH", "HHHBPH", "HHHHBB", "HHPBHH",
                                 "HPBHHP", "PHHBPH","trail")
) {
    read_tsv(mffile,comment="#",skip=1,
             col_names=mitofates_colnames) %>%
        dplyr::select(-trail) %>%
        mutate(Pred_preseq=factor(Pred_preseq,
                                  levels=c("No mitochondrial presequence",
                                           "Possessing mitochondrial presequence"),
                                  labels=c("No","Yes")))  %>%
        mutate(Gene = str_remove(Gene, " Protein"))
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

infonfromseqs <- function(seqs,sstart=1L,send=-1L,ssize=4,pseudocount=0.5,
                          scope=c("tot","all")) {
    Cmat <- seqs %>% 
    str_sub(start=sstart,end=send) %>%
    str_to_upper %>%
    consensusMatrix
    if (scope == "tot") {
        Cmat %>%
            infontot(ssize=ssize,pseudocount=pseudocount) %>%
            return()
    } else     if (scope == "all") {
        Cmat %>%
            infonall(ssize=ssize,pseudocount=pseudocount) %>%
            return()
    } else {
        warning("invalid scope option to infonfromseqs")
    }
}

plotheat_inf_perbase <- function(inf_perbase,
                                         Keylevels=c("All","HiTrans","CytoRibo"),
                                 inf_limits=c(0,2)) {
    ggplot(data=inf_perbase %>% 
               gather(key=Key,value="Infon",-i,-Pos) %>%
               mutate(Key=factor(Key, levels=Keylevels )),
           aes(x=i,y=Key,fill=Infon) ) +
        geom_tile() + 
        scale_fill_gradient(low="white",high="darkblue",limits=inf_limits) + 
        theme(axis.title=element_blank())
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

PWMscoren <- function(seqs,pwm=kozak_n_PWM,startl=8) {
    PWMscore(seqs,pwm,startl) 
}

PWMscorew <- function(seqs,pwm=kozak_w_PWM,startl=4) {
    PWMscore(seqs,pwm,startl) 
}

is.context <- function(x) str_detect(x,"context")


