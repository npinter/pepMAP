# pepMAP

pepMAP is a web-based application for visualizing peptide mappings onto protein sequences. It supports input from both FragPipe and DIA-NN reports, providing an interactive interface to explore peptide coverage, modifications, and protein features.

![pepMAP](https://github.com/user-attachments/assets/ddbb7964-f953-48a8-b094-c35d5491d49d)

## Features

- **Support for FragPipe and DIA-NN Inputs**: Upload peptide reports from FragPipe or DIA-NN to visualize peptide mappings on proteins.

- **Interactive Visualization**: Explore peptide coverage with an intuitive interface powered by Plotly.

- **Protein Feature Integration**: Fetch and display protein domains, binding sites, modified residues, and variants from UniProt.

- **Multi-User Sessions**: Supports multiple users simultaneously with isolated sessions.

- **Custom Features Upload**: Accepts custom features either with the initial upload or in a later step.

## Installation

To run pepMAP locally, follow these steps:

### 1. Clone the Repository

```bash
git clone https://github.com/npinter/pepMAP.git
cd pepMAP
```

### 2. Create a Conda Environment

Create a new Conda environment named `pepMAP` with the required packages:

```bash
conda create -n pepMAP python=3.9
conda activate pepMAP
```

### 3. Install Dependencies

Install the necessary packages with the specified versions:

```bash
conda install -c conda-forge numpy==1.26.4 pandas==2.2.2 plotly==5.21.0 requests==2.31.0 flask==3.0.3 flask-caching==2.1.0 apscheduler==3.10.4
```

### 4. Run the Application

Activate the Conda Environment

```bash
conda activate pepMAP
```
Start the Flask application:

```bash
(pepMAP) python app.py
```

By default, the application runs on `http://localhost:7007`. Open this URL in your web browser to access pepMAP.

## Usage

1. **Upload Files**: Upload your peptide report (from FragPipe or DIA-NN) and FASTA file used in the search.

2. **Search Proteins**: Enter a UniProt ID or gene symbol to visualize peptide mappings.

3. **Custom Features Format**: Use a tab-separated file with the columns `uniprot`, `position`, `description`, `literature`. The `position` field can be a single number or a range.

Example (`custom_features.tsv`):

```tsv
uniprot position    description literature
P04626  1242-1254   IHC (Ventana 4B5)   https://doi.org/10.1111/j.1365-2559.2011.04034.x
P04626  557-603 ADC (Trastuzumab/Herceptin) https://doi.org/10.1038/nature01392
Q9NZQ7  284-290 IHC (Ventana SP263) https://doi.org/10.1038/s41379-019-0372-z
```


## License

This project is licensed under the MIT License.


# Note

- Ensure that you have an active internet connection, as the application fetches protein features from the UniProt API.
- The application stores uploaded data in `storage/` and schedules daily deletion to manage server resources effectively.
