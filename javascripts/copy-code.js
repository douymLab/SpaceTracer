document.addEventListener("DOMContentLoaded", function () {
  const pres = document.querySelectorAll(".rst-content pre");

  pres.forEach((pre) => {
    if (pre.dataset.copyButtonAdded === "true") return;
    pre.dataset.copyButtonAdded = "true";

    const code = pre.querySelector("code");
    if (!code) return;

    const wrapper = document.createElement("div");
    wrapper.className = "copy-btn-container";

    pre.parentNode.insertBefore(wrapper, pre);
    wrapper.appendChild(pre);

    const button = document.createElement("button");
    button.className = "copy-code-button";
    button.type = "button";
    button.setAttribute("aria-label", "Copy code to clipboard");
    button.innerHTML = `
      <svg viewBox="0 0 24 24" aria-hidden="true">
        <path d="M16 1H4c-1.1 0-2 .9-2 2v12h2V3h12V1zm3 4H8c-1.1 0-2 .9-2 2v14c0 1.1.9 2 2 2h11c1.1 0 2-.9 2-2V7c0-1.1-.9-2-2-2zm0 16H8V7h11v14z"></path>
      </svg>
    `;

    button.addEventListener("click", async () => {
      try {
        await navigator.clipboard.writeText(code.innerText);
        const originalHTML = button.innerHTML;
        // Replace with check icon
        button.innerHTML = `
          <svg viewBox="0 0 24 24" aria-hidden="true">
            <path d="M9 16.2l-3.5-3.5L4 14.2l5 5 12-12-1.5-1.5z"></path>
          </svg>
        `;
        button.classList.add("copied");

        setTimeout(() => {
          button.innerHTML = originalHTML;
          button.classList.remove("copied");
        }, 1200);
      } catch (err) {
      console.error("Copy failed:", err);
      }
    });

    wrapper.appendChild(button);
  });
});
