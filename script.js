document.addEventListener("DOMContentLoaded", function () {
    const sections = document.querySelectorAll(".body-section");
    const navLinks = document.querySelectorAll(".sidebar .content-section");
    const rightSidebarContent = document.getElementById("right-sidebar-content");

    const sidebarMessages = {
        "intro-section": "Additional information for introduction section. Some vocab and clarifications for what SNPs, genotypes, phenotypes, etc are.",
        "data-section": "Links to data(?).",
        "methods-section": "Additional information about our methods. Definitions for LD and Thresholding, or we can put the flowchart here.",
        "results-section": "Any additional notes for the results section.",
        "additional-content-section": "Find extra resources and discussions here.",
    };

    // Navigate to each section by clicking on the sidebar with smooth scrolling
    navLinks.forEach(anchor => {
        anchor.addEventListener("click", function (e) {
            e.preventDefault(); // Prevent default link behavior

            const targetId = this.getAttribute("href").substring(1); // Remove #
            const targetElement = document.getElementById(targetId);

            if (targetElement) {
                window.scrollTo({
                    top: targetElement.offsetTop - 60, // Adjust scroll position
                    behavior: "smooth"
                });
            }
        });
    });

    function highlightCurrentSection() {
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

        // Update right sidebar content with fade effect
        if (sidebarMessages[currentSection]) {
            rightSidebarContent.style.opacity = "0"; // Fade out
            setTimeout(() => {
                rightSidebarContent.textContent = sidebarMessages[currentSection]; // Change text
                rightSidebarContent.style.opacity = "1"; // Fade in
            }, 300);
        }
    }

    window.addEventListener("scroll", highlightCurrentSection);
});
