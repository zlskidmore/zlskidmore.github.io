document.addEventListener('DOMContentLoaded', () => {
  const toggleSwitch = document.getElementById('theme-toggle'); // the checkbox input
  const root = document.documentElement;

  const storedTheme = localStorage.getItem('theme') || 'dark';
  root.setAttribute('data-theme', storedTheme);
  toggleSwitch.checked = storedTheme === 'light';

  // Toggle logic
  toggleSwitch.addEventListener('change', () => {
    const newTheme = toggleSwitch.checked ? 'light' : 'dark';
    root.setAttribute('data-theme', newTheme);
    localStorage.setItem('theme', newTheme);
  });
});