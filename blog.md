---
layout: default
title: Blog
permalink: /blog/
pagination:
  enabled: true
  per_page: 1
---

<div class="container">
  <header class="mb-4">
    <h1>Blog</h1>
    <p class="lead">Latest posts and notes.</p>
  </header>

  <div class="list-group">
    {% for post in paginator.posts %}
      <article class="list-group-item">
        <h3 class="mb-1"><a href="{{ post.url | relative_url }}">{{ post.title }}</a></h3>
        <div class="small text-muted mb-2">
          {{ post.date | date: "%B %d, %Y" }} · {% capture words %}{{ post.content | number_of_words }}{% endcapture %}{% assign minutes = words | plus:0 | divided_by:200 %}
          {% if minutes < 1 %}☕ Quick read{% else %}⏱ {{ minutes }} min read{% endif %}
          {% if post.tags %} · Tags:
            {% for tag in post.tags %}
              <a href="{{ '/tag/' | append: tag | slugify | prepend: '/' | replace: '//','/' }}">{{ tag }}</a>{% unless forloop.last %}, {% endunless %}
            {% endfor %}
          {% endif %}
        </div>
        <p class="mb-0">{{ post.excerpt | strip_html | truncatewords: 40 }}</p>
      </article>
    {% endfor %}
  </div>

  <!-- Pagination controls -->
<nav aria-label="Blog pagination" class="my-4">
  <ul class="pagination justify-content-center">
    {% if paginator.previous_page %}
      <li class="page-item">
        <a class="page-link" href="{{ paginator.previous_page_path | relative_url }}">← Newer</a>
      </li>
    {% else %}
      <li class="page-item disabled"><span class="page-link">← Newer</span></li>
    {% endif %}

    {% for page_num in (1..paginator.total_pages) %}
      {% assign path = paginator.paginate_path | replace: ':num', page_num %}
      {% if page_num == paginator.page %}
        <li class="page-item active"><span class="page-link">{{ page_num }}</span></li>
      {% else %}
        <li class="page-item"><a class="page-link" href="{{ path | relative_url }}">{{ page_num }}</a></li>
      {% endif %}
    {% endfor %}

    {% if paginator.next_page %}
      <li class="page-item">
        <a class="page-link" href="{{ paginator.next_page_path | relative_url }}">Older →</a>
      </li>
    {% else %}
      <li class="page-item disabled"><span class="page-link">Older →</span></li>
    {% endif %}
  </ul>
</nav>
</div>
