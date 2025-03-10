document.addEventListener("DOMContentLoaded", function () {
    const sections = document.querySelectorAll(".body-section");
    const navLinks = document.querySelectorAll(".sidebar .content-section");
    const rightSidebarContent = document.getElementById("right-sidebar-content");

    let hoveringPlot = false; // Flag to track hover state

    const sidebarMessages = {
        "intro-section": "Single Nucleotide Polymorphism (SNP): genetic variant at a specific location \n\n Genotype: unique DNA code; a set of instructions that tells your body what to do \n\n Gene expression: the set of instructions actually being used. Can change over time due to factors like environment, age, and disease. \n\n Genome-wide association studies (GWAS): identify genetic variants associated with a disease \n\n Transcriptome-wide association studies (TWAS): identify genes associated with a disease via their gene expression",
        "data-section": "LDREF: Linkage Disequilibrium Reference",
        "methods-section": "Linkage Disequilibrium (LD): non-random association of alleles at two or more loci \n\n Thresholding: setting a threshold to determine significant SNPs",
        "results-section": "",
        "additional-content-section": "",
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

    addHoverEffect(document.querySelector(".loci-container"),
        "Loci of significant genes:<br><br>We used the loci of significant genes to assess genetic similarities in heart failure suspectibility. As demonstrated on the plot, the amount of overlap between the populations is pretty modest, suggesting that the genes associated with heart disease are very population-specific."
    );

    addHoverEffect(document.querySelector(".huge-container"),
        "HuGE ORMDL3 Gene Association:<br><br>Using a stricter significance threshold of p < 0.001, we identified several genes strongly associated with heart failure in the European population. We determined the strength of the associations using HuGE (Human Genetic Evidence) Score. Specifically, gene ORMDL3 in chromosome 17 exhibited the strongest association with heart failure, as it had the highest HuGE score for this condition. "
    );

    addHoverEffect(document.querySelector(".go-container"), 
        "Gene Ontology (GO) Enrichment Analysis:<br><br>GO Enrichment Analysis is a method used to find which biological functions are most affected by the significant genes in each population. For instance, we can infer from the plot that hypoxia (low oxygen) is a major factor in heart failure in the European population, followed by regulation of smooth muscle contraction. <br><br>On the other hand, the biological processes displayed for the American population do not exhibit that strong of an association with heart failure comparing to the European population."
    );
});

function toggleDropdown(dropdownId, button) {
    let dropdown = document.getElementById(dropdownId);

    if (dropdown.style.display === "block") {
        dropdown.style.display = "none";
        button.textContent = "Click for more info";
    } else {
        dropdown.style.display = "block";
        button.textContent = "Click to hide";
    }
}




