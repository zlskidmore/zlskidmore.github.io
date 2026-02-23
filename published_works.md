---
layout: default
title: Published Works
permalink: /published_works/
---

### Overview of My Research and Publications

{% assign total = 0 %}
{% for year in site.data.publications %}
  {% assign total = total | plus: year[1].size %}
{% endfor %}

Welcome to my publications page.

You'll find a strong focus on cancer genomics, spanning everything from individual case reports to large-scale genomic studies and clinical trials. A few papers also highlight my contributions to bioinformatics software development. As of now, I have {{ total }} published papers — and there's always more to come. I aim to keep this list current, so feel free to check back periodically.

– Zach

A complete list of publications can be found on:
<a href="https://www.ncbi.nlm.nih.gov/myncbi/zach.skidmore.1/bibliography/public/">
PubMed <i class="fa-solid fa-book"></i>
</a>
 | 
<a href="https://www.ncbi.nlm.nih.gov/myncbi/zach.skidmore.1/bibliography/public/">
ORCiD <i class="fa-brands fa-orcid"></i>
</a>
 | 
<a href="https://www.ncbi.nlm.nih.gov/myncbi/zach.skidmore.1/bibliography/public/">
Google Scholar <i class="fa-solid fa-file-lines"></i>
</a>

<hr>

<div class="container my-5">
  <div class="accordion" id="pubAccordion">
    {% include accordion_pub_year.html year="2026" papers=site.data.publications.y2026 %}
    {% include accordion_pub_year.html year="2025" papers=site.data.publications.y2025 %}
    {% include accordion_pub_year.html year="2024" papers=site.data.publications.y2024 %}
    {% include accordion_pub_year.html year="2023" papers=site.data.publications.y2023 %}
    {% include accordion_pub_year.html year="2022" papers=site.data.publications.y2022 %}
    {% include accordion_pub_year.html year="2021" papers=site.data.publications.y2021 %}
    {% include accordion_pub_year.html year="2020" papers=site.data.publications.y2020 %}
    {% include accordion_pub_year.html year="2019" papers=site.data.publications.y2019 %}
    {% include accordion_pub_year.html year="2018" papers=site.data.publications.y2018 %}
    {% include accordion_pub_year.html year="2017" papers=site.data.publications.y2017 %}
    {% include accordion_pub_year.html year="2016" papers=site.data.publications.y2016 %}
    {% include accordion_pub_year.html year="2015" papers=site.data.publications.y2015 %}
  </div>
</div>