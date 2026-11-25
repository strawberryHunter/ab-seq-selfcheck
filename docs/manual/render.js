// Take screenshots for the user manual and print it to PDF.
//
// This file is for building the docs. abcheck itself uses no external libraries;
// this script (Playwright) is needed only when rebuilding the manual PDF.
//
//   node docs/manual/render.js
//
// See docs/manual/build.sh for the order of steps.

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

  // 1) Report screenshots -- used in chapter 6 of the manual
  const shotPage = await browser.newPage({
    viewport: { width: 1000, height: 1400 },
    deviceScaleFactor: 2,
  });
  await shotPage.goto(REPORT);
  await shot(shotPage, ".summary", "report_summary.png");
  await shot(shotPage, "#seq-SYN-VH-004", "report_detail.png");
  await shot(shotPage, "#seq-SYN-VH-901", "report_failed.png");
  await shotPage.close();

  // 2) Manual PDF
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
