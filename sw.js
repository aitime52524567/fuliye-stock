/* ============================================================
 * A股周期分析 PWA Service Worker
 * 作用：缓存应用外壳，实现「加到主屏幕后」可离线启动
 * 数据接口（东方财富行情/搜索）不走缓存，始终联网
 * ============================================================ */

const CACHE = 'astock-cycle-v1';
const ASSETS = [
  './',
  './index.html',
  './style.css',
  './app.js',
  './echarts.min.js',
  './manifest.json',
  './icon-192.png',
  './icon-512.png',
  './icon-180.png'
];

self.addEventListener('install', (e) => {
  e.waitUntil(
    caches.open(CACHE)
      .then((c) => c.addAll(ASSETS))
      .then(() => self.skipWaiting())
  );
});

self.addEventListener('activate', (e) => {
  e.waitUntil(
    caches.keys()
      .then((keys) => Promise.all(
        keys.filter((k) => k !== CACHE).map((k) => caches.delete(k))
      ))
      .then(() => self.clients.claim())
  );
});

self.addEventListener('fetch', (e) => {
  if (e.request.method !== 'GET') return;
  const url = new URL(e.request.url);

  // 行情 / 搜索接口：网络直连，不缓存（保证数据实时）
  if (url.hostname.endsWith('eastmoney.com')) return;

  // 其余资源（含 echarts CDN）：缓存优先 + 后台更新
  e.respondWith(
    caches.match(e.request).then((cached) => {
      const network = fetch(e.request).then((res) => {
        if (res && (res.status === 200 || res.type === 'opaque')) {
          const clone = res.clone();
          caches.open(CACHE).then((c) => c.put(e.request, clone));
        }
        return res;
      }).catch(() => cached);
      return cached || network;
    })
  );
});
