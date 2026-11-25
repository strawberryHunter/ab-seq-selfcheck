// 사용 설명서의 스크린샷을 찍고 PDF 로 인쇄한다.
//
// 이 파일은 문서 빌드용이다. abcheck 자체는 외부 라이브러리를 쓰지 않으며,
// 이 스크립트(Playwright)는 설명서 PDF 를 다시 만들 때만 필요하다.
//
//   node docs/manual/render.js
//
// 순서는 docs/manual/build.sh 참고.

const path = require("path");
const { chromium } = require("playwright");

const ROOT = path.resolve(__dirname, "..", "..");
const REPORT = "file://" + path.join(ROOT, "sample_output", "report.html");
const MANUAL = "file://" + path.join(__dirname, "manual.html");
const IMG = path.join(__dirname, "img");
const PDF = path.join(__dirname, "abcheck_manual.pdf");

async function shot(page, selector, file) {
  const el = page.locator(selector).first();
  await el.scrollIntoViewIfNeeded();
  await el.screenshot({ path: path.join(IMG, file) });
  console.log("스크린샷:", file);
}

(async () => {
  const browser = await chromium.launch();

  // 1) 리포트 스크린샷 — 매뉴얼 6장에 들어간다
  const shotPage = await browser.newPage({
    viewport: { width: 1000, height: 1400 },
    deviceScaleFactor: 2,
  });
  await shotPage.goto(REPORT);
  await shot(shotPage, ".summary", "report_summary.png");
  await shot(shotPage, "#seq-SYN-VH-004", "report_detail.png");
  await shot(shotPage, "#seq-SYN-VH-901", "report_failed.png");
  await shotPage.close();

  // 2) 매뉴얼 PDF
  const page = await browser.newPage();
  await page.goto(MANUAL, { waitUntil: "load" });
  await page.pdf({
    path: PDF,
    format: "A4",
    printBackground: true,
    preferCSSPageSize: true,
    displayHeaderFooter: true,
    headerTemplate: "<span></span>",
    footerTemplate:
      '<div style="width:100%;font-size:8px;color:#777;padding:0 17mm;' +
      'display:flex;justify-content:space-between;font-family:sans-serif;">' +
      "<span>abcheck 사용 설명서</span>" +
      '<span><span class="pageNumber"></span> / <span class="totalPages"></span></span></div>',
  });
  console.log("PDF:", path.relative(ROOT, PDF));
  await browser.close();
})();
