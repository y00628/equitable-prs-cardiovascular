document.addEventListener("DOMContentLoaded", function () {
    const sections = document.querySelectorAll("section");
    const navLinks = document.querySelectorAll(".content-section");

    // Smooth scrolling
    navLinks.forEach(anchor => {
        anchor.addEventListener("click", function (e) {
            e.preventDefault();

            const targetId = this.getAttribute("href").substring(1);
            const targetElement = document.getElementById(targetId);

            if (targetElement) {
                window.scrollTo({
                    top: targetElement.offsetTop - 20, // Adjusts scroll position
                    behavior: "smooth"
                });
            }
        });
    });

    // Active section highlighting
    function updateActiveSection() {
        let scrollPosition = window.scrollY;

        sections.forEach(section => {
            const sectionTop = section.offsetTop - 50; // Adjust for header height
            const sectionHeight = section.offsetHeight;
            const sectionId = section.getAttribute("id");

            if (scrollPosition >= sectionTop && scrollPosition < sectionTop + sectionHeight) {
                navLinks.forEach(link => {
                    link.classList.remove("active");
                });
                document.querySelector(`.content-section[href="#${sectionId}"]`).classList.add("active");
            }
        });
    }

    window.addEventListener("scroll", updateActiveSection);
    updateActiveSection(); // Call on page load to highlight the correct section
});
