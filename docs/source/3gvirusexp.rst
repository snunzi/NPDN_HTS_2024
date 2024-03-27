In-Person 3G Viral Metagenomics
===================================


Introduction
^^^^^^^^^^^^
Here, we will analyze the data you generated over the past two days. You should each have two samples, one from a tomato host and the other from an orange host. Lets see what viruses we have!


Import Data
^^^^^^^^^^^
Lets import data from a shared history. These are your raw reads. The only thing that has been done is the multiple read files have been concatenated into a single file for analysis. There is also a database of plant viruses, and your two host genomes.

.. admonition:: Hands-On: Import Viral Metagenomic Reads

    1. At the top of the screen click on ``Shared Data`` and select ``Histories``

    2. In the search field, search for ``NPDN HTS2 2024``

    3. Find the history for ``NPDN HTS2 2024 Data 3G Virus `` Select the button at the top to import into your Galaxy environment.

    4. After uploading the shared history you can delete barcodes that are not yours and rename your sequences to identify what they are, i.e. 'tomato_metagenome.fastq.gz' and 'orange_metagenome.fastq.gz'.

    5. You should now have 2 sequencing files (fastq.gz), two host genomes (orange and tomato), and a virus database (nt_vrl_plant_pathoscope_95.fasta) in your history.



Sequence QC
^^^^^^^^^^^^^
The first step in any sequencing analysis is quality check and trimming. These sequences have already been based called with on-board base calling and this is how you would receive them off of the sequencer. Let's first check the quality of the data we received.


.. admonition:: Hands-On: Quality Check

    1. In tools menu, search for 'Nanoplot' and click on it.

    2. Run Nanoplot tool with the following parameters

    * files: Select "Multiple Datasets" optin

    * “files”: ``tomato_metagenome.fastq.gz`` and ``citrus_metagenome.fastq.gz``

    * Leave the rest as default.

    3. Click Run Tool.


Nanoplot should produce four output files. Let's take a look at the html output report. Fill in your handout with QC information and select your read length cutoff for filtering.


Quality Filtering
^^^^^^^^^^^^^^^^^^^
Let's filter the data to remove short reads and low quality bases. First we will filter to retain only high-quality long reads. Quality filtering is a balancing act to retain enough high-quality reads for analysis. Here, we will set a minimum length for reads to maintain. We will also only keep the top 90% of high quality reads.



.. admonition:: Hands-On: Quality Filtering

    1. In tools menu, search for 'Filtlong' and click on it.

    2. Run Filtlong tool with the following parameters

      * Input Fastq: ``porechop output``

      * Output Theshholds:

          - Min Mean Quality: ``10``

          - Min Length: ``1000`` (Or cutoff you selected)

      * Leave the rest as default.

    3. Click Execute.

    4. After filtering completes, run Nanoplot again to see how many reads remain and their distribution. Fill in handout with output.





Non-Host Read Extraction
^^^^^^^^^^^^^^^^^^^^^^^^^^

We will now remove host reads from the dataset. We will use the two host genomes you imported from the shared history. These are the RefSeq genomes from tomato and orange. You will have to run the following two programs (minimap and samtools) separately for each sample because they are from different host genomes.

.. admonition:: Hands-On: Remove Host Reads

    1. Run minimap2 with the following parameters:

      * Will you select a reference genome from your history or use a built-in index? ``Use a genome from history and build index``

      * Use the following dataset as the reference sequence: choose either citrus or tomato (based on which sample)

      * Select fastq datasets: ``filtlong output``

      * Leave rest as default press 'Run Tool'


Now we have to extract only the reads that do not map to the host genome.
    2. Run samtools fastx

      * “BAM or SAM file to convert”: ``Map with minimap2``

      * “Output format”: ``compressed FASTQ``

      * “Outputs”: ``others``

      * “Require that these flags are set”: ``Read is unmapped``

      * Leave rest as default press 'Run Tool'

    3. When job completes, rename the output files to something more useful.

      * Click on pencil icon next to ``data X converted to fastqsanger.gz`` and rename to ``meta_tomato_nonhost.fastq.gz``


Read Assignment with minimap2
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We will be mapping all reads to identify members in a mixed set of metagenomic reads.

.. admonition:: Hands-On: Viral Read Mapping with minimap2


    1. Run minimap2 with the following parameters:

      * Will you select a reference genome from your history or use a built-in index?: ``Use a genome from history and build index``

      * Use the following dataset as the reference sequence:  ``nt_vrl_plant_pathoscope_95.fasta``

      * Select fastq dataset: ``tomato_nonhost.fastq.gz`` and ``citrus_nonhost.fastq.gz``

    2. Run tool.

.. admonition:: Hands-On: Count read Mapping

    1. Find tool ``samtools idxstats``

    2. Run samtools idxstats with the following parameters:

    * BAM file: ``Map with minimap...`` (you can select both files--one for tomato and one for citrus)

    3. Run tool.

    4. Download this file to your computer and open in excel to examine.

Extra-Pull out mapped Reads
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. admonition:: Hands-On: Pull out mapped reads

    1. Make a bed file of what you need, this is a text file with name of your genome (tab) 1	(tab)	length of your genome:

  * Formated:

  * ti|1761477|gi|0|ref|OM515245|Tomato_brown_rugose_fruit_virus_isolate_39986372,_complete_genome	1	7767

  2. Import your bed file into Galaxy

  3. Open the tool ``samtools view``. Run with the following parameters:

  * SAM/BAM/CRAM data set : ``Minimap bam file``

  * What would you like to look at? ``A filtered/subsampled selection of reads``

  * Configure Filters:

  * Filter by regions: ``Regions from a BED file``

  * Filter by intervals in a bed file ``Select your bed file``

  * Run Tool

  4. Open the tool ``samtools fastx``. Run with the following parameters:

    * “BAM or SAM file to convert”: ``Filtered bam``

    * “Output format”: ``fasta``

  5. You can blast these reads to a larger database (NCBI) to see what they are



Metagenome Assembly
^^^^^^^^^^^^^^^^^^^^^

Next we will assemble all reads that did not map to host using an assembler for 3G data, Flye. There are multiple assemblers available for MinION data, but this assembler provides a nice balance of accuracy and speed.

.. admonition:: Hands-On: Assembly with Flye

    1. In the tools menu search for 'flye' tool and click on it.

    2. Run this tool with following parameters:

      * Input Reads: ``citrus_nonhost.fastq.gz`` and ``tomato_nonhost.fastq.gz``


      * Perform metagenomic assembly: ``Yes``

      * Leave the rest as default

    3. Run tool.

When the assembly completes, take a look at the ``Flye assembly info`` output.




Blast Contigs
^^^^^^^^^^^^^^

Let's Blast the contigs we generated. First we will build a blast database

.. admonition:: Hands-On: Blast Contigs

  1. Find the tool ``NCBI BLAST+ makeblastdb``.

  2. Run this tool with the following parameters:

  * Molecule type of input: ``nucleotide``

  * Subject database/sequences: ``Blast database from your history``

  3. Run tool.



.. admonition:: Hands-On: Blast Contigs

  1. Find the tool ``NCBI BLAST+ blastn``.

  2. Run this tool with the following parameters:

  * Nucleotide query sequence(s): ``Flye Consensus``

  * Subject database/sequences: ``Blast database from your history``

  3. Run tool.

  4. Download results to computer and open in excel.

  Coverage of Assemblies
  ^^^^^^^^^^^^^^^^^^^^^^^

  Let's look at the coverage of the assemblies.

  .. admonition:: Hands-On: Turn Genome Assembly into a Database

    1. Follow this tutorial to turn your genome assemblies into databases: https://training.galaxyproject.org/training-material/faqs/galaxy/analysis_add_custom_build.html

    2. Next to Database in the assembly history panel, click the question mark and assign the database.

      .. image:: _static/databasegalaxy.png

    3. Search for the tool minimap2 and run with the following parameters;

    * Will you select a reference genome from your history or use a built-in index?: ``Use a genome from history and build index``

    * Use the following dataset as the reference sequence:  ``tomato assembly`` (or citrus)

    * Select fastq dataset: ``tomato_nonhost.fastq.gz`` (or``citrus_nonhost.fastq.gz``)

    3. Run tool.

    4. When mapper completes select 'visualize' in the history panel.
