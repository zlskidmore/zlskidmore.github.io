document.addEventListener("DOMContentLoaded", function () {
  const posts = document.querySelectorAll('.post-card');
  const tagFilterContainer = document.getElementById('tag-filters');
  const allTags = new Set();

  posts.forEach(post => {
    const tags = post.dataset.tags.split(' ');
    tags.forEach(tag => allTags.add(tag));
  });

  // Create "All" button
  const allBtn = document.createElement('button');
  allBtn.className = 'btn btn-outline-dark m-1 active';
  allBtn.textContent = 'All';
  allBtn.dataset.filter = '*';
  tagFilterContainer.appendChild(allBtn);

  // Create tag buttons
  Array.from(allTags).sort().forEach(tag => {
    const btn = document.createElement('button');
    btn.className = 'btn btn-outline-dark m-1';
    btn.textContent = tag;
    btn.dataset.filter = tag;
    tagFilterContainer.appendChild(btn);
  });

  // Add filter behavior
  tagFilterContainer.addEventListener('click', function (e) {
    if (e.target.tagName !== 'BUTTON') return;

    const filter = e.target.dataset.filter;

    document.querySelectorAll('#tag-filters .btn').forEach(btn => btn.classList.remove('active'));
    e.target.classList.add('active');

    posts.forEach(post => {
      const postTags = post.dataset.tags.split(' ');
      const show = filter === '*' || postTags.includes(filter);
      post.style.display = show ? 'block' : 'none';
    });
  });
});
