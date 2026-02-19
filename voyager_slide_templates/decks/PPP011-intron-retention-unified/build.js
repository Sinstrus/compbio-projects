const { html2pptx } = require('../../html2pptx/index.cjs');
const PptxGenJS = require('pptxgenjs');
const fs = require('fs');
const path = require('path');

async function buildPresentation() {
  console.log('Building Unified Intron Retention VHH Display presentation...\n');

  const pres = new PptxGenJS();
  pres.layout = 'LAYOUT_16x9';
  pres.title = 'Intron Retention VHH Display on AAV9';
  pres.author = 'TRACER-Nano';
  pres.subject = 'Unified Intron Retention: VR4 Library + N-Terminal Positions';
  pres.company = 'Voyager Therapeutics';

  const slides = [
    "slide01-title.html",
    "slide02-executive-summary.html",
    "slide03-dual-product.html",
    "slide04-cassette-architecture.html",
    "slide05-cryptic-sites.html",
    "slide06-u1-framework.html",
    "slide07-tiered-donors.html",
    "slide08-gene-derived.html",
    "slide09-vr4-verification.html",
    "slide10-why-nterm.html",
    "slide11-avd132-vp2.html",
    "slide12-avd133-vp1.html",
    "slide13-splicing-topology.html",
    "slide14-risk-verification.html",
    "slide15-conclusions.html"
  ];

  for (let i = 0; i < slides.length; i++) {
    const slideNum = String(i + 1).padStart(2, '0');
    const slidePath = path.join(__dirname, slides[i]);

    console.log(`  Processing slide ${slideNum}...`);

    try {
      await html2pptx(slidePath, pres, {
        tmpDir: fs.realpathSync(path.join(__dirname, 'assets')),
      });
      console.log(`    ✓ Slide ${slideNum} added`);
    } catch (error) {
      console.error(`    ✗ Error on slide ${slideNum}:`, error.message);
      throw error;
    }
  }

  console.log(`\nTotal slides: ${slides.length}`);
  console.log('Saving presentation...\n');

  try {
    const outputPath = path.join(__dirname, 'Intron_Retention_VHH_Display_Unified.pptx');
    await pres.writeFile({ fileName: outputPath });

    const stats = fs.statSync(outputPath);
    console.log(`SUCCESS: Presentation saved to ${outputPath}`);
    console.log(`File size: ${(stats.size / 1024).toFixed(1)} KB`);
  } catch (error) {
    console.error('ERROR saving file:', error.message);
    process.exit(1);
  }
}

buildPresentation();
