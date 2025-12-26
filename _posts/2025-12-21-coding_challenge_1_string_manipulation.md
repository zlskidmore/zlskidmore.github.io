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

Solutions in both R and Python follow each challenge. They’re not the only way to solve the problem, just clean, reproducible examples to compare your own work against.

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
"TCGACCGTT"
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
find all k-mers

### Challenge 5
count base frequencies

### Challenge 6
codon grouping

### Challenge 8
Compute GC content in 10bp sliding windows

### Challenge 9
count base mismatches between forward and reverse dna strands of the same molecule

### Challenge 10
Find the motif

