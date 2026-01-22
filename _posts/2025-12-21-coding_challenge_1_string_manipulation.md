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

Solutions in both R and Python follow each challenge. They’re not the only way to solve the problem, just clean, reproducible examples to compare your own work against. They are in fact the solutions I came up with on the fly, so you might notice a bias in libraries, using data.table in R for example :)

Have Fun!

<hr>

### Challenge 1 - Validate DNA Sequence

DNA sequences should consist only of the four canonical bases: A, C, G, and T. In practice, errors or data corruption can introduce invalid characters.

Given a DNA string, identify any invalid bases and report their 1-based positions along with the incorrect characters.

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

### Challenge 2 - Reverse Complement a DNA Sequence

DNA is conventionally written in the 5′→3′ direction, but many analyses require working with the opposite strand.

Given a DNA sequence, compute its reverse complement, ensuring the output is returned in uppercase.

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
comp = str.maketrans({"A": "T", "T": "A", "C": "G", "G": "C"})

rev_comp = dnaSeq.upper().translate(comp)[::-1]
print(rev_comp)
```
</details>

<hr>

### Challenge 3 - Calculate GC Content

The relative abundance of G and C bases can have a real biological influence due to the three hydrogen bonds in GC base pairs, making GC-rich regions more thermodynamically stable.

Given a DNA string, calculate the GC content as a percentage of the total sequence length.

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
seq = dnaSeq.upper()
gc_pct = 100 * sum(1 for b in seq if b in ("G", "C")) / len(seq)
print(gc_pct)
```
</details>

<hr>

### Challenge 4 - Generate All Possible k-mers

k-mers, short sequences of length k, are a core concept in genomics, appearing in alignment, assembly, and indexing algorithms.

For a given value of k, generate all possible DNA k-mers composed of A, C, G, and T, and report the total number of unique combinations.

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
from itertools import product

kmers = ["".join(p) for p in product(possibleBases, repeat=k)]
print(len(kmers))
```
</details>

<hr>

### Challenge 5 - Base Composition

Understanding the relative composition of nucleotide bases is a common first step in exploratory sequence analysis.

Given a DNA string, calculate the proportion of each base (A, C, G, T) within the sequence.

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
seq = dnaSeq.lower()
counts = pd.Series(list(seq)).value_counts()
props = (counts / counts.sum()).reindex(["a", "c", "g", "t"]).fillna(0)

out = pd.DataFrame({"bases": props.index, "proportions": props.values})
print(out)
```
</details>

<hr>


### Challenge 6 - Open Reading Frames

Protein-coding regions begin at a start codon (ATG) and terminate at the first encountered stop codon (TAA, TAG, or TGA).

Given a DNA sequence, locate the first ATG codon and split the sequence into codons (groups of three bases) from that point until the first stop codon is reached.

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
# setup
dnaSeq = "ACTTTCTTATGTTTAGTTTCAATATTGTTTTCTTTTCTCTGGCTAATAAAGGCCTTATTCATTTCTAATTATGAAA"
```

###### Solution

<details markdown="1">
  <summary>Show</summary>

```python
import re

m = re.search("ATG", dnaSeq)
if not m:
    codons = []
else:
    orf = seq[m.start():]
    codons = [orf[i:i+3] for i in range(0, len(orf), 3)]

    stops = {"TAA", "TAG", "TGA"}
    stop_i = next((i for i, c in enumerate(codons) if c in stops), None)
    if stop_i is not None:
        codons = codons[:stop_i + 1]

print(codons)
```
</details>

<hr>

### Challenge 7 - Finly poly-a tracts

Homopolymer runs, such as poly-A tracts, are common features of biological sequence data and can represent true biological signals (e.g., poly-A tails) or technical artifacts.

Given a DNA sequence written 5′→3′, identify all poly-A runs consisting of three or more consecutive A’s. For each run, report the 1-based start position, end position, and run length. Runs should be maximal, meaning each run should extend as far as possible. Ignore lowercase versus uppercase characters.

###### Output

```R
   start   end    sequence
   <int> <int>      <char>
1:    19    24      AAAAAA
2:    30    40 AAAAAAAAAAA
3:    52    54         AAA
```

#### R

###### Input

```R
# setup
dnaSeq <- "TCGTGCCTGACGCAATGCAAAAAAGTCGCAAAAAAAAAAATGGCTGCGCTCAAA"
```
###### Solution

<details markdown="1">
  <summary>Show</summary>

```R
# load lib
library(stringr)
library(data.table)

# find positions of all poly-a runs, 3 or more A's in a row
poly_a_pos <- as.data.table(str_locate_all(dnaSeq, "A{3,}"))

# with positions extracted we can just build the sequence back
poly_a_pos[,sequence := str_dup("A", (end-start) + 1)]
```
</details>

<br>

#### Python

###### Input

```python
# setup
dnaSeq = "TCGTGCCTGACGCAATGCAAAAAAGTCGCAAAAAAAAAAATGGCTGCGCTCAAA"
```

###### Solution

<details markdown="1">
  <summary>Show</summary>

```python
import re
import pandas as pd

runs = []
for m in re.finditer(r"A{3,}", dnaSeq):
    start = m.start() + 1   # 1-based
    end = m.end()           # end is already 1-based if we treat end as inclusive
    runs.append((start, end, m.group()))

df = pd.DataFrame(runs, columns=["start", "end", "sequence"])
print(df)
```
</details>

<hr>

### Challenge 8 - GC Content in a Sliding Window

GC content often varies across a sequence, and examining it locally can reveal regions of low complexity or unusual composition.

Given a DNA string, compute the GC content across a 10-base sliding window. The window should be right-aligned, such that the first window evaluated corresponds to bases 1–10, followed by bases 2–11, and so on.

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
# setup
dnaSeq = "ACTTTCTTATGTTTAGTTTCAATATTGTTTTCTTTTCTCTGGCT"
```

###### Solution

<details markdown="1">
  <summary>Show</summary>

```python
import pandas as pd

gc01 = [1 if b in ("G", "C") else 0 for b in dnaSeq]

window = 10
s = pd.Series(gc01)

# rolling defaults to right-aligned; require full window like your NA filtering
gc_pct = (s.rolling(window=window).mean() * 100).dropna()

# print as ints if you want to match your output style closely
print(gc_pct.astype(int).tolist())
```
</details>

<hr>

### Challenge 9 - Base Pair Mismatch

Base mismatches between complementary DNA strands may represent true variants or technical artifacts.

Given two DNA strings, each written 5′→3′, reverse one strand to align it antiparallel to the other. Then compare bases using standard complement rules (A - T, G - C) and count the number of mismatches.

For example:

5'-ATGCC-3'
3'-TACGT-5'

Contains one artifact at the very end on the right side, either the C or T is wrong.

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
# setup
strand_1 = "ATGCCGTCA"
strand_2 = "ACACTGCAT"
```

###### Solution

<details markdown="1">
  <summary>Show</summary>

```python
comp = {"A": "T", "T": "A", "C": "G", "G": "C"}

s1 = strand_1.upper()
s2 = strand_2.upper()[::-1]  # reverse to align antiparallel

mismatches = sum(1 for a, b in zip(s1, s2) if comp.get(a) != b)
print(mismatches)
```
</details>

<hr>

### Challenge 10 - Locate a Motif

Short sequence motifs play important roles in gene regulation and genome annotation.

Given a DNA sequence and a target motif, we'll uses a TATA-box motif, identify the 1-based start and end positions of the motif within the sequence.

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
dnaSeq = "CGCTATAAAAGGGC"
tataBox = "TATAAAAG"
```

###### Solution

<details markdown="1">
  <summary>Show</summary>

```python
import re
import pandas as pd

m = re.search(re.escape(tataBox), dnaSeq)

if not m:
    df = pd.DataFrame(columns=["start", "end"])
else:
    start = m.start() + 1
    end = m.end()          # inclusive end in 1-based terms
    df = pd.DataFrame([(start, end)], columns=["start", "end"])

print(df)
```
</details>

<hr>
