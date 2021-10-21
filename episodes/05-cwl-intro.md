[Go to main page](../README.md)

## Intro to workflows.


### What is Common Workflow Language (CWL) ?

CWL is a way to describe command line tools and connect them together to create workflows. Because CWL is a specification and not a specific piece of software, tools and workflows described using CWL are portable across a variety of platforms that support the CWL standard.

### Why might I want to learn to use CWL ?

CWL has roots in “make” and many similar tools that determine order of execution based on dependencies between tasks. However unlike “make”, CWL tasks are isolated and you must be explicit about your inputs and outputs. The benefit of explicitness and isolation are flexibility, portability, and scalability: tools and workflows described with CWL can transparently leverage technologies such as Docker and be used with CWL implementations from different vendors. CWL is well suited for describing large-scale workflows in cluster, cloud and high performance computing environments where tasks are scheduled in parallel across many nodes.

### First example

The simplest "hello world" program.  This accepts one input parameter, writes a message to the terminal or job log, and produces no permanent output. CWL documents are written in [JSON][json] or [YAML][yaml], or a mix of the two. We will use YAML throughout this guide. If you are not familiar with YAML, you may find it helpful to refer to [this quick tutorial for the subset of YAML used in CWL](http://www.commonwl.org/user_guide/yaml/).

First, create a file called `1st-example.cwl`, containing the boxed text below. It will help you to use a text editor that can be specified to produce text in YAML or JSON. Whatever text editor you use, the indents you see should not be created using tabs.

*1st-tool.cwl*

~~~
{% include cwl/1st-example.cwl %}
~~~
{: .source}

Next, create a file called `1st-example.yml`, containing the following boxed text, which will describe the input of a run:

*echo-job.yml*

~~~
{% include cwl/1st-example.yml %}
~~~
{: .source}

Now, invoke `cwl-runner` (or `cwltool`) with the tool wrapper `1st-example.cwl` and the input object echo-job.yml on the command line. The command
is  `cwltool 1st-example.cwl 1st-example.yml`. The boxed text below shows this command and the expected output.

~~~
$ cwl-runner 1st-example.cwl 1st-example.yml
[job 1st-tool.cwl] /tmp/tmpmM5S_1$ echo \
    'Hello world!'
Hello world!
[job 1st-tool.cwl] completed success
{}
Final process status is success
~~~
{: .output}