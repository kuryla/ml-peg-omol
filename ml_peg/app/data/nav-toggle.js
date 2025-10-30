// Simple mobile nav toggle for Dash app header
document.addEventListener('DOMContentLoaded', function () {
  const btn = document.getElementById('nav-toggle');
  const nav = document.getElementById('site-nav');
  if (!btn || !nav) return;
  btn.addEventListener('click', function () {
    const isOpen = nav.classList.toggle('open');
    btn.setAttribute('aria-expanded', isOpen ? 'true' : 'false');
  });
  // highlight active link based on location
  const links = nav.querySelectorAll('a.nav-link');
  const path = window.location.pathname;
  links.forEach(a => {
    a.classList.remove('active');
    if ((path === '/' && a.getAttribute('href') === '/') ||
        (path.startsWith('/models') && a.getAttribute('href') === '/models') ||
        (path.startsWith('/datasets') && a.getAttribute('href') === '/datasets')) {
      a.classList.add('active');
    }
  });
});
