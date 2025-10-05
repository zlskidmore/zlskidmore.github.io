const toggleBtn = document.getElementById('theme-toggle');
const root = document.documentElement;

// Set default theme to dark unless stored otherwise
const storedTheme = localStorage.getItem('theme') || 'dark';
root.setAttribute('data-theme', storedTheme);

toggleBtn.addEventListener('click', () => {
  const currentTheme = root.getAttribute('data-theme');
  const newTheme = currentTheme === 'dark' ? 'light' : 'dark';
  root.setAttribute('data-theme', newTheme);
  localStorage.setItem('theme', newTheme);
});

