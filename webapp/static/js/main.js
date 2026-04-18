// Main JS - shared utilities
// Currently all page-specific JS is inline in templates
// This file reserved for future shared utilities

document.addEventListener('DOMContentLoaded', () => {
  // Mark current nav link as active
  const path = window.location.pathname;
  document.querySelectorAll('.nav-links a').forEach(a => {
    if (a.getAttribute('href') === path) {
      a.classList.add('active');
    }
  });
});
