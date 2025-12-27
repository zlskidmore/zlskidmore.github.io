---
layout: post
title:  "Coding Challenges #1: String Manipulation with DNA Sequences"
date:   2025-12-21 10:13:50 -0500
description: "Ten mini coding exercises involving sequence manipulation in R and Python"
tags: [Code_Challenge, R, Python]
---

### Overview

Welcome to the first post in a new series of bite-sized coding challenges.

Each challenge is intentionally small: just enough to keep your skills sharp without being annoying. The goal is simple, stay fluent in R and Python by solving practical, self-contained problems.

For every prompt, I'll provide clearly defined input and a target output. Sometimes you may read data from a file; other times you’ll generate it on the fly. Either way, the idea is the same: run the setup code, then try to reproduce the expected result.

Solutions in both R and Python follow each challenge. They’re not the only way to solve the problem, just clean, reproducible examples to compare your own work against. They are in fact the solutions I came up with on the fly, so you might notice a bias twoards using data.table in R for example :)

Have Fun!

<hr>

### Challenge 1

DNA base sequences consist of A, C, G, and T. Find the wrong bases and print out their positions!

###### Output

```R
"Indices: 9, 14 contain incorrect bases: H, Z"
```

#### R

###### Setup

```R
# Setup
dnaSeq <- "ACGTCCTGHGCACZTTA"
```

###### Solution

<details markdown="1">
  <summary>Show</summary>

```R
# split string into a vector, find the indices, and then extract the bases
dnaVec <- unlist(strsplit(dnaSeq, ""))
badBase_idx <- which(!dnaVec %in% c('A', 'C', 'G', 'T'))
badBase <- dnaVec[badBase_idx]

# format the output
print(paste0("1-Based Indices: ", toString(badBase_idx), " contain incorrect bases: ", toString(badBase)))
```
</details>

<br>

#### Python

###### Input

```python
# Setup
dnaSeq = "ACGTCCTGHGCACZTTA"
```

###### Solution

<details markdown="1">
  <summary>Show</summary>

```python
# define a "set" of good bases
goodBases = {"A", "C", "G", "T"}

# use string enumeration, read like this, [WHAT_TO_KEEP for ITEM in ITERABLE if CONDITION]
# WHAT_TO_KEEP is: i
# ITEM is: i, base
# ITERABLE is: enumerate(dnaSeq)
indexes = [i for i, base in enumerate(dnaSeq) if base not in goodBases]

# Run through iterable again, use the index to extract bases
bad_bases = [dnaSeq[i] for i in indexes]

# format for printing
print(f"Indices: {indexes} contain incorrect bases: {bad_bases}")
```
</details>

<hr>

### Challenge 2

DNA is by convention written 5' - 3', a somewhat intuitive exercise is to display the opposite strand of a dna fragment, in other words reverse complement the dna strand. Make sure everything is upper case.

###### Output

```R
"TCGACCGTT"
```

#### R

###### Input

```R
# Setup
dnaSeq <- "AACggtCGA"
```
###### Solution

<details markdown="1">
  <summary>Show</summary>

```R
# convert to uppercase
dnaSeq_rc <- toupper(dnaSeq)

# use switch statement to change bases
switchBase <- function(base) {
  
  base <- as.character(base)
  compBase <- switch(base, "A"="T", "T"="A", "C"="G", "G"="C")
  return(compBase)
}
seq_vec <- unlist(strsplit(dnaSeq_rc, ""))
seq_vec <- sapply(seq_vec, switchBase)

# reverse sequence
dnaSeq_rc <- paste0(rev(seq_vec), sep="", collapse="")
dnaSeq_rc
```
</details>

<br>

#### Python

###### Input

```python
# Setup
dnaSeq = "AACggtCGA"
```

###### Solution

<details markdown="1">
  <summary>Show</summary>

```python

```
</details>

<hr>

### Challenge 3

The ratio of G and C to A and T in DNA sequences can introduce bias, calculate the GC content of a string of bases as a percentage.

###### Output

```R
48.14815
```

#### R

###### Input

```R
# Setup
dnaSeq <- "acgcgtcgacgttttgccataatatcg"
```
###### Solution

<details markdown="1">
  <summary>Show</summary>

```R
# load lib
library(data.table)

# split up into a character vector of bases
dnaSeq_vec <- unlist(strsplit(dnaSeq, ""))

# convert to a data.table
dnaSeq_dt <- as.data.table(dnaSeq_vec)
setnames(dnaSeq_dt, "bases")

# count the overal base counts and calculate GC content
dnaSeq_counts <- dnaSeq_dt[,.N,by=.(bases)]
answer <- sum(dnaSeq_counts[bases %in% c('g', 'c')]$N)/sum(dnaSeq_counts$N) * 100
```
</details>

<br>

#### Python

###### Input

```python
# Setup
dnaSeq = "acgcgtcgacgttttgccataatatcg"
```

###### Solution

<details markdown="1">
  <summary>Show</summary>

```python

```
</details>

<hr>

### Challenge 4

Occasionally you may be tasked with finding all possible DNA sequences for a specific k-mer length, try that here, for a kmer of length 6 (i.e. 6 bases in length), find all possible DNA sequence combinations.
find all k-mers

###### Output

```R
4096
```

#### R

###### Input

```R
# setup
possibleBases <- c("A","C","G","T")
k <- 6
```
###### Solution

<details markdown="1">
  <summary>Show</summary>

```R
# we can use expand.grid, which will generate all possible combinations for a list of vectors
# we create a list here using just a simply lambda function
a <- function(x) {
  return(possibleBases)
}
basesList <- lapply(1:k, a)

# calling expand.grid with the list produces all combination
grid <- expand.grid(basesList, stringsAsFactors = FALSE)

# anonymous function to collapse columns for each row
a <- function(x){
  paste0(x, sep="", collapse="")
}
kmer_list <- apply(grid, 1, a)
length(kmer_list)
```
</details>

<br>

#### Python

###### Input

```python
# setup
possibleBases <- ['A', 'C', 'G', 'T']
k = 6
```

###### Solution

<details markdown="1">
  <summary>Show</summary>

```python

```
</details>

<hr>

### Challenge 5

From a DNA string, count the proportion of each base within the sequence, don't worry about printing the exact output, just get the answer for a, c, g, and t

###### Output

```R
    bases    proportions
1:      a      0.2222222
2:      c      0.2592593
3:      g      0.2222222
4:      t      0.2962963
```

#### R

###### Input

```R
# setup
dnaSeq <- "acgcgtcgacgttttgccataatatcg"
```
###### Solution

<details markdown="1">
  <summary>Show</summary>

```R
# load lib
library(data.table)

# split the string
dnaVec <- unlist(strsplit(dnaSeq, ""))

# convert to a data.table
dnaDT <- as.data.table(dnaVec)
setnames(dnaDT, 'bases')

# construct counts and proportions
dnaCountDT <- dnaDT[,.N, by=.(bases)]
dnaCountDT[,baseProportion := N/sum(N)]
```
</details>

<br>

#### Python

###### Input

```python
# Setup
dnaSeq = "acgcgtcgacgttttgccataatatcg"
```

###### Solution

<details markdown="1">
  <summary>Show</summary>

```python

```
</details>

<hr>


### Challenge 6

Translation of DNA begins at methionine amino acid encoded in DNA by the letters ATG, find the ATG codon (groups of 3 DNA bases) and identify the remaining codons up to the stop site (TAA, TAG, or TGA).

###### Output

```R
"ATG" "TTT" "AGT" "TTC" "AAT" "ATT" "GTT" "TTC" "TTT" "TCT" "CTG" "GCT" "AAT" "AAA" "GGC" "CTT" "ATT" "CAT" "TTC" "TAA"
```

#### R

###### Input

```R
# setup
dnaSeq <- "ACTTTCTTATGTTTAGTTTCAATATTGTTTTCTTTTCTCTGGCTAATAAAGGCCTTATTCATTTCTAATTATGAAA"
```
###### Solution

<details markdown="1">
  <summary>Show</summary>

```R
# load lib
library(stringr)
library(data.table)

# find the start and strip out whats before, use a regex and non-greedy quantifier with a capture group,
# this represents the start of the ORF
dnaSeq_orf <- gsub("^.*?(ATG.*)", "\\1", dnaSeq)

# split into codons using str_sub
start <- seq(1, nchar(dnaSeq_orf), by=3)
stop  <- pmin(start + 2, nchar(dnaSeq_orf))
dna_Seq_codons <- str_sub(dnaSeq_orf, start=start, end=stop)

# annotate stop codons
dna_Seq_codons_DT <- data.table(codon=dna_Seq_codons)
dna_Seq_codons_DT[,StopCodon := grepl("^(TAA|TAG|TGA)$", codon)]

# subset to pull up until the first stop encounter, then format back to a vector
dna_Seq_codons_DT <- dna_Seq_codons_DT[1:min(which(dna_Seq_codons_DT$StopCodon == TRUE))]
dna_Seq_codons_DT$codon
```
</details>

<br>

#### Python

###### Input

```python

```

###### Solution

<details markdown="1">
  <summary>Show</summary>

```python

```
</details>

<hr>

### Challenge 8

In this challenge we will examine the GC content of a string of DNA across a rolling window, given the DNA string below, compute the GC content for a 10bp sliding window. The window should be right aligned so the first complete window to calculate on would be "ACTTTCTTAT" then "CTTTCTTATG" etc.

###### Output

```R
20 30 20 20 20 20 20 20 20 20 30 20 20 20 20 20 10 20 20 20 10 10 20 20 20 20 20 20 20 30 30 40 40 50 50
```

#### R

###### Input

```R
# setup
dnaSeq <- "ACTTTCTTATGTTTAGTTTCAATATTGTTTTCTTTTCTCTGGCT"
```
###### Solution

<details markdown="1">
  <summary>Show</summary>

```R
# lib
library(data.table)

# plan is to use DT frollapply for a rolling window, we need an integer vector for that, so
# we will re-encode the data, G and C == 1 and A and T == 0
dnaSeq_vec <- unlist(strsplit(dnaSeq, ""))
dnaSeq_vec <- ifelse(grepl("G|C", dnaSeq_vec), 1, 0)

# create a function for GC content calculation as a percentage
a <- function(x) {
  gc_content <- sum(x)/length(x) * 100
  return(gc_content)
}

# apply over a 10 bp window, right align the window is the default
windowedGC <- frollapply(dnaSeq_vec, 10, a)
windowedGC[!is.na(windowedGC)]
```
</details>

<br>

#### Python

###### Input

```python

```

###### Solution

<details markdown="1">
  <summary>Show</summary>

```python

```
</details>

<hr>

### Challenge 9

Commonly, base mismatches occur in bioinformatics data. These may represent true variants or technical artifacts. Below are two DNA strings, each written 5′→3′. To compare them, reverse one strand to align it antiparallel to the other, then compare bases using standard complement rules and count the number of mismatches.
For example:

5'-ATGCC-3'
3'-TACGT-5'

Contains one artifact at the very end on the right side, either the C or T is wrong.

Remember:
G matches C
A matches T

###### Output

```R
3
```

#### R

###### Input

```R
strand_1 <- "ATGCCGTCA"
strand_2 <- "ACACTGCAT"
```
###### Solution

<details markdown="1">
  <summary>Show</summary>

```R
# lib
library(data.table)

# make a DT to hold the forward and reverse DNA strands
dnaSeq_DT <- data.table(forward = unlist(strsplit(strand_1, "")),
                        reverse = rev(unlist(strsplit(strand_2, ""))))

# examine each outcome and annotate if it's a bad alignment
dnaSeq_DT[forward == 'A', badAlign := ifelse(reverse != "T", 1, 0)]
dnaSeq_DT[forward == 'T', badAlign := ifelse(reverse != "A", 1, 0)]
dnaSeq_DT[forward == 'C', badAlign := ifelse(reverse != "G", 1, 0)]
dnaSeq_DT[forward == 'G', badAlign := ifelse(reverse != "C", 1, 0)]

# count the number of bad alignments
sum(dnaSeq_DT$badAlign)
```
</details>

<br>

#### Python

###### Input

```python

```

###### Solution

<details markdown="1">
  <summary>Show</summary>

```python

```
</details>

<hr>

### Challenge 10

Often in biology we are interested in motif's, below I supply a DNA sequence, find the 1-base start and stop position of the TATA box for the sequence in the setup

###### Output

```R
     start end
[1,]     4  11
```

#### R

###### Input

```R
dnaSeq <- "CGCTATAAAAGGGC"
tataBox <- "TATAAAAG"
```
###### Solution

<details markdown="1">
  <summary>Show</summary>

```R
# lib
library(stringr)

# locate the pattern
str_locate(dnaSeq, tataBox)
```
</details>

<br>

#### Python

###### Input

```python

```

###### Solution

<details markdown="1">
  <summary>Show</summary>

```python

```
</details>

<hr>
