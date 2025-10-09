document.addEventListener("DOMContentLoaded", function () {
  const posts = document.querySelectorAll('.post-card');
  const tagFilterContainer = document.getElementById('tag-filters');
  const tagCounts = {};

  // Build tag count map
  posts.forEach(post => {
    const tags = post.dataset.tags.split(' ');
    tags.forEach(tag => {
      tagCounts[tag] = (tagCounts[tag] || 0) + 1;
    });
  });

  // Create "All" button
  const allBtn = document.createElement('button');
  allBtn.className = 'btn btn-outline-dark m-1 active';
  allBtn.textContent = `All (${posts.length})`;
  allBtn.dataset.filter = '*';
  tagFilterContainer.appendChild(allBtn);

  // Create one button per tag, sorted alphabetically
  Object.keys(tagCounts).sort().forEach(tag => {
    const btn = document.createElement('button');
    btn.className = 'btn btn-outline-dark m-1';
    btn.textContent = `${tag} (${tagCounts[tag]})`;
    btn.dataset.filter = tag;
    tagFilterContainer.appendChild(btn);
  });

  // Core filtering logic
  tagFilterContainer.addEventListener('click', function (e) {
    if (e.target.tagName !== 'BUTTON') return;

    const filter = e.target.dataset.filter;

    // Toggle active button style
    document.querySelectorAll('#tag-filters .btn').forEach(btn => {
      btn.classList.remove('active');
    });
    e.target.classList.add('active');

    posts.forEach(post => {
      const postTags = post.dataset.tags.split(' ');
      const shouldShow = filter === '*' || postTags.includes(filter);

      // First: if card is hidden but should be shown, prepare it
      if (shouldShow && post.classList.contains('hidden')) {
        post.classList.remove('hidden', 'is-hiding');
        // Force reflow to ensure animation runs
        void post.offsetWidth;
        post.style.display = 'block'; // just in case
      }

      // Now apply transitions
      if (!shouldShow && !post.classList.contains('is-hiding')) {
        post.classList.add('is-hiding');
        // After animation ends, fully hide the element
        setTimeout(() => {
          post.classList.add('hidden');
        }, 1000); // Match CSS duration
      }
    });
  });
});
