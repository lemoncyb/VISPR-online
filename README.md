# VISPR-online
VISPR-online is a web-based interactive framework for CRISPR screens visualization, exploration and sharing.

#### Table of Contents
1. [Description](#Description)
2. [Installation](#Installation)
3. [Demo Test](#Demo)
4. [FAQs](#FAQs)
5. [License](#License)
6. [Contact](#Contact)

<a name="Description"></a>
I. Description
----
VISPR-online is a web-based interactive framework for CRISPR screens visualization, exploration and sharing.

CRISPR screening helps systematically exploration of the functions of coding and non-coding elements in a genome. We previously developed [MAGeCK](https://www.ncbi.nlm.nih.gov/pubmed/25476604) and [MAGeCK-VISPR](https://www.ncbi.nlm.nih.gov/pubmed/26673418) to perform CRISPR screening data analysis. VISPR is a visualization tool included in [MAGeCK-VISPR](https://www.ncbi.nlm.nih.gov/pubmed/26673418), which can be used to explore interesting genes. However, VISPR is designed for local use and some manual configutations need be set to run the program. To help the community to access the tool easily, we improve the tool and develop an online version: VISPR-online.

The advantages of VISPR-online compared with VISPR:

* Installation and configuration free. Easy to access via a browser.
* Support interactive view of sgRNA locations in the gene context.
* Enable to resume and share data and sessions.

<a name="Installation"></a>
II. Installation
----
VISPR-online can also be installed in a local computer or network for internal use.

### Step1:
Download VISPR-online source code.

```
$ git clone https://github.com/lemoncyb/VISPR-online.git
```

### Step2:
To install VISPR-online you have to use the Python 3 variant of the Miniconda Python distribution (http://conda.pydata.org/miniconda.html). VISPR-online cannot be installed on the Python 2 variant of Miniconda. When installing Miniconda, make sure that you answer yes to this question:
```
Do you wish the installer to prepend the Miniconda3 install location to PATH ...? [yes|no]
```
Also, make sure that you do not have set the PYTHONPATH environment variable, because it will interfere with the Miniconda setup.

### Step3:
Enter the top directory of VISPR-online:
```
$cd VISPR-online
```
On the terminal or an Anaconda Prompt, you can create an isolated software environment for vispr-online by executing:
```
$conda env create -f ./environment.yml
```
***Note:***
If the error *ResolvePackageNotFound:-pandas==0.19.1* is reported, please skip to *[FAQs](#FAQs)* setup.

### Step4:
Then, activate the environment by running
```
$conda activate vispr-online
```
conda activate and conda deactivate only work on conda 4.6 and later versions. For conda versions prior to 4.6, run:

Windows: 
```
$activate vispr-online
```

Linux and macOS: 
```
$source activate vispr-online
``` 

  For detailed manual of miniconda, please refer to https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html 


<a name="Demo"></a>
III. Demo test
----
### Step1:
Launch VISPR-online server.

Enter the top directory of VISPR-online. Execute the following script:

Windows:
```
$python ./run.py
```
Linux and macOS: 
```
$ ./run.py
```

Open your browser and access <http://127.0.0.1:5000>. You will see the following page if successful.

![](vispr_screen/static/img/homepage.jpg) 

### Step2:
A demo is integrated in VISPR-online repository. You could click the demo link on home page to load the project and explore the results.

### Step3:
A test dataset is contained in VISPR-online repository. The path is ```VISPR-online/vispr_screen/static/data/testdata.zip```. The dataset contains the following files:

* **mle.gene_summary.txt:** gene summary file.
* **all.count_normalized.txt:** normalized count file.
* **mle.sgrna_summary.txt:** sgRNA summary file.
* **sgrnas.bed:** sgRNA location file.
* **README.txt:** README file.

***Note:***
`The input of VISPR-online are direct output of screening analysis tools, so the input file type and file name are customized for different tools. User could match files by the file suffix. For example, the name of gene summary file of MAGeCK ends with “.gene_summary.txt”. Besides, the corresponding file suffix will appear when the mouse hovers over the file type, as shown in the figure below.`

![](vispr_screen/static/img/file_type.jpg) 

We summarize the file suffix of VISPR-online input in the table below. We have added this table and more detailed information of input file in the tutorial page.
![](vispr_screen/static/img/file_type_table.jpg) 

Upload corresponding files to VISPR-online. Then click the "Submit" button to explore the results. The results pages are like this.

![](vispr_screen/static/img/result_view.jpg)  

<a name="FAQs"></a>
IV. FAQs
----
### Error1: 
If the error *ResolvePackageNotFound:-pandas==0.19.1* is reported, please follow the steps below to install:

#Step1: Create an environment with a specific version of Python:
```
$conda create -n vispr-online python=3.5
```
#Step2: Activate the environment by running
```
$conda activate vispr-online
```
#Step3: Install related Python dependencies.
```
$ pip install flask
$ pip install pymongo
$ pip install PyYAML
$ pip install numpy
$ pip install pandas==0.19.1 -i http://pypi.douban.com/simple --trusted-host pypi.douban.com
$ pip install sklearn
```
### Error2:
If the error *ImportError: cannot import name 'EmptyDataError' from 'pandas.io.common'* is reported, it is because the version of panda is not compatible. Please follow the above steps to install pandas *version 0.19.1*.

### Error3:
The error *"/usr/bin/env: 'python': No such file or directory"* when running the run.py script. This is because the first line of the run.py script assumes that python is an alias to python3. Please install python3 to resolve this error.

<a name="License"></a>
License
----
Licensed under the [MIT license](http://opensource.org/licenses/MIT). This project may not be copied, modified, or distributed except according to those terms.

<a name="Contact"></a>
Contact
----
Yingbo Cui <yingbocui@nudt.edu.cn>
