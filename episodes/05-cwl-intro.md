[Go to main page](../README.md)

## Intro to workflows.


### What is Common Workflow Language (CWL) ?

CWL is a way to describe command line tools and connect them together to create workflows. Because CWL is a specification and not a specific piece of software, tools and workflows described using CWL are portable across a variety of platforms that support the CWL standard.

### Why might I want to learn to use CWL ?

CWL has roots in “make” and many similar tools that determine order of execution based on dependencies between tasks. However unlike “make”, CWL tasks are isolated and you must be explicit about your inputs and outputs. The benefit of explicitness and isolation are flexibility, portability, and scalability: tools and workflows described with CWL can transparently leverage technologies such as Docker, be used with CWL implementations from different vendors, and is well suited for describing large-scale workflows in cluster, cloud and high performance computing environments where tasks are scheduled in parallel across many nodes.

### Overview of CWL syntax 

~~~
$ cwl-runner 1st-tool.cwl echo-job.yml
[job 1st-tool.cwl] /tmp/tmpmM5S_1$ echo \
    'Hello world!'
Hello world!
[job 1st-tool.cwl] completed success
{}
Final process status is success

~~~