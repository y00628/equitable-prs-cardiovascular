document.addEventListener("DOMContentLoaded", function () {
    const sections = document.querySelectorAll(".body-section");
    const navLinks = document.querySelectorAll(".sidebar .content-section");
    const rightSidebarContent = document.getElementById("right-sidebar-content");

    let hoveringPlot = false; // Flag to track hover state

    const sidebarMessages = {
        "intro-section": "Single Nucleotide Polymorphism (SNP): genetic variant at a specific location \n\n Genotype: unique DNA code; a set of instructions that tells your body what to do \n\n Gene expression: the set of instructions actually being used. Can change over time due to factors like environment, age, and disease. \n\n Genome-wide association studies (GWAS): identify genetic variants associated with a disease \n\n Transcriptome-wide association studies (TWAS): identify genes associated with a disease via their gene expression",
        "data-section": "Links to data(?).",
        "methods-section": "Linkage Disequilibrium (LD): non-random association of alleles at two or more loci \n\n Thresholding: setting a threshold to determine significant SNPs",
        "results-section": "Hover over a plot to view more insights",
        "additional-content-section": "Find extra resources and discussions here.",
    };

    // Function to update the right sidebar text
    function updateSidebarContent(content) {
        if (!hoveringPlot) { // Only update if no plot is being hovered
            rightSidebarContent.style.opacity = "0"; // Fade out
            setTimeout(() => {
                rightSidebarContent.innerHTML = content.replace(/\n/g, "<br>"); // Preserve new lines
                rightSidebarContent.style.opacity = "1"; // Fade in
            }, 200);
        }
    }

    // Smooth scrolling for sidebar navigation
    navLinks.forEach(anchor => {
        anchor.addEventListener("click", function (e) {
            e.preventDefault();
            const targetId = this.getAttribute("href").substring(1);
            const targetElement = document.getElementById(targetId);

            if (targetElement) {
                window.scrollTo({
                    top: targetElement.offsetTop - 60,
                    behavior: "smooth"
                });
            }
        });
    });

    // Highlight the current section in the sidebar
    function highlightCurrentSection() {
        if (hoveringPlot) return; // Skip if a plot is being hovered

        let maxVisibleHeight = 0;
        let currentSection = "";

        sections.forEach((section) => {
            const rect = section.getBoundingClientRect();
            const visibleHeight = Math.min(window.innerHeight, rect.bottom) - Math.max(0, rect.top);

            if (visibleHeight > maxVisibleHeight) {
                maxVisibleHeight = visibleHeight;
                currentSection = section.getAttribute("id");
            }
        });

        // Update sidebar link highlighting
        navLinks.forEach((link) => {
            link.classList.remove("active");
            if (link.getAttribute("href").includes(currentSection)) {
                link.classList.add("active");
            }
        });

        // Update right sidebar content
        if (sidebarMessages[currentSection]) {
            updateSidebarContent(sidebarMessages[currentSection]);
        }
    }

    window.addEventListener("scroll", highlightCurrentSection);

    // Handle hover behavior for plots
    function addHoverEffect(container, text) {
        container.addEventListener("mouseover", function () {
            hoveringPlot = true; // Set hover state
            rightSidebarContent.style.opacity = "0"; // Fade out
            setTimeout(() => {
                rightSidebarContent.innerHTML = text;
                rightSidebarContent.style.opacity = "1"; // Fade in
            }, 200);
        });

        container.addEventListener("mouseout", function () {
            hoveringPlot = false; // Reset hover state
            highlightCurrentSection(); // Restore the default section-based sidebar text
        });
    }

    // Define hover content for each plot
    addHoverEffect(document.querySelector(".miami-container"),
        "Miami Plot of GWAS Results:<br><br>This Miami plot presents genome-wide association study (GWAS) results, comparing two populations: African ancestry (AFR, top half) vs. European ancestry (EUR, bottom half).<br><br>Some SNPs appear significant in one population but not the other, suggesting genetic variants that may influence cardiovascular disease risk differently between African and European ancestries. For instance, LINC00565 is significant in the European population, but not the African population."
    );

    addHoverEffect(document.querySelector(".manhattan-container"),
        "EAS Manhattan Plot:<br><br>The Manhattan plot visualizes the statistical significance of genetic variants associated with a trait across the genome. Peaks indicate regions with significant associations. A few SNPs exceed the genome-wide significance threshold (dashed line), suggesting potentially important genetic associations in this population. Notably, significant SNPs appear on chromosomes 4, 9, 10, and 12, indicating candidate regions for further investigation.<br><br>"
    );

    addHoverEffect(document.querySelector(".venn-container"),
        "Venn Diagram of Population Gene Overlaps:<br><br>Using the Venn diagram, we compared the significant genes in each population using the European cohort as a reference, given that SNP weights were initially trained on this population. We found out that only 3 to 4 overlapping genes per population. The modest degree of overlap reinforces the importance of population-specific PRS models."
    );

    addHoverEffect(document.querySelector(".huge-container"),
        "Huge ORMDL3 Gene Association:<br><br>ORMDL3 has been strongly linked to inflammatory responses and respiratory diseases. This visualization explores its association with disease susceptibility."
    );
});
