## 📚 I. Onboarding & Community Guidelines

This section sets the tone and provides initial context.

* **Code of Conduct:** Clearly state expectations for respectful and professional communication (often links to a separate `CODE_OF_CONDUCT.md`).
* **Getting Started:**
    * **Prerequisites:** List required software (Python 3.x, Docker, Git).
    * **Installation:** Provide the exact steps to clone the repository and run the application locally using `docker-compose up` or the local development setup (mentioning the virtual environment, e.g., `env`).
* **Support:** Where to ask questions (e.g., GitHub Discussions, a dedicated Slack/Discord channel).

## 🛠️ II. Plugin Development Technical Specifications

This is the most critical section for your architecture. It defines the "contract" for a valid plugin.

### A. Plugin Structure and Registration
* **The Template:** Direct contributors to copy a simple, pre-built template directory (e.g., `omni_bio_ai/plugin_template`) to start their work.
* **Required Files:** Every plugin **must** contain:
    * `models.py`: For database structures (if any).
    * `views.py`/`urls.py`: For its web endpoints.
    * `plugin_config.py`: Must inherit from `omni_bio_ai.core.BaseBioPlugin` and define mandatory attributes:
        * `plugin_name` (e.g., `'crispr_analysis'`)
        * `display_title` (e.g., `'CRISPR Off-Target Predictor'`)
        * `icon_class` (e.g., `'fa-dna'`)
* **Registration:** Explain how the plugin is registered (e.g., "Add the plugin's name to the `PLUGINS` list in the main `settings.py` or, ideally, register using **Python Entry Points**").

### B. Data & Inter-Plugin Communication
* **The Central Data Model:**
    * Plugins **must not** store raw data locally. They must interact with the **Central Data Lake/Broker** (`data_uploader` app).
    * **Input/Output Format:** Define standard data transfer formats (e.g., all genomic data is accessed via a **File ID**; all simple output should be structured JSON or a standard bioinformatics format like VCF/BAM/FASTQ).
* **Calling Other Plugins:** Specify the standard for inter-plugin API calls (e.g., "Use the `omni_bio_ai.services.api_facade` to call another plugin's REST endpoint").

### C. LLM/AI Integration
* **LLM Interface:** Define the canonical way to request AI inference (e.g., "All LLM calls must go through the `ollama-server` service via the `rag_inference` API, ensuring proper model tracking and resource limits are respected").
* **Custom ML Models:** Explain where researchers should place custom ML model files (e.g., in the `models_store` directory) and how they should be loaded by the `ml_predictor` component.

## ✍️ III. Code Quality and Testing

Since you plan to scale, strict quality control is essential.

* **Coding Style:** Enforce a standard (e.g., **PEP 8** for Python, mention use of code formatters like **Black** or **isort**).
* **Testing Requirement (Mandatory):**
    * **Unit Tests:** Every plugin **must** include its own self-contained test suite in its dedicated `tests/` directory.
    * **Coverage:** Require a minimum test coverage percentage (e.g., 80% or 90%) for the added code.
    * **Testing Mocks:** Provide clear examples of how to **mock** external services (like the LLM, database, or API calls) to ensure unit tests are fast and isolated.
* **Documentation:** All public functions, classes, and models must have **Docstrings** following a specific format (e.g., NumPy or Google style).

## 🚀 IV. Submission Process

This defines how a contribution moves from development to live code.

1.  **Issue Creation:** **Must** create a GitHub issue first to discuss the new plugin idea and get approval on scope.
2.  **Branching:** Use a clear branch naming convention (e.g., `feature/plugin-name` or `fix/issue-number`).
3.  **Pull Request (PR) Checklist:** Include a non-negotiable checklist for the PR description:
    * [ ] Passes local unit tests.
    * [ ] Includes documentation for the new feature.
    * [ ] Addresses the issue [Link to Issue #].
    * [ ] Provides a test dataset for validation (if applicable).
    * [ ] Commits are signed (if you use a DCO/CLA).
4.  **Review and Merge:** Set expectations on review time and the need for reviewers to verify the plugin's integration with the core platform.

These guidelines ensure that new plugins are not only functional but are architecturally sound, testable, and maintainable within the larger OmniBioAI ecosystem.
