[Go to main page](../README.md)

## How do I run tools or workflows in parallel?

Now that we know how to write workflows, we can start utilizing the `ScatterFeatureRequirement`. This feature tells the runner that you wish to run a tool or workflow multiple times over a list of inputs. The workflow then takes the input(s) as an array and will run the specified step(s) on each element of the array as if it were a single input. This allows you to run the same workflow on multiple inputs without having to generate many different commands or input yaml files.

~~~
requirements:
  ScatterFeatureRequirement: {}
~~~
{: .source}

The most common reason a new user might want to use scatter is to perform the same analysis on
different samples. Let's start with a simple workflow that calls our first example and takes
an array of strings as input to the workflow:

***scatter-workflow.cwl***

~~~
{% include cwl/scatter-workflow.cwl %}
~~~
{: .source}

Aside from the `requirements` section including `ScatterFeatureRequirement`, what is
going on here?

~~~
inputs:
  message_array: string[]
~~~
{: .source}

First of all, notice that the main workflow level input here requires an array of strings.

~~~
steps:
  echo:
    run: echo.cwl
    scatter: message
    in:
      message: message_array
    out: []
~~~
{: .source}

Here we've added a new field to the step `echo` called `scatter`. This field tells the
runner that we'd like to scatter over this input for this particular step. Note that
the input name listed after scatter is the one of the step's input, not a workflow level input.

For our first scatter, it's as simple as that! Since our tool doesn't collect any outputs, we
still use `outputs: []` in our workflow, but if you expect that the final output of your
workflow will now have multiple outputs to collect, be sure to update that to an array type
as well!

Using the following input file:

***scatter.yml***

~~~
{% include cwl/scatter-job.yml %}
~~~
{: .source}

As a reminder, `1st-tool.cwl` simply calls the command `echo` on a message. If we invoke
`cwl-runner scatter-workflow.cwl scatter-job.yml` on the command line:

~~~
$ cwl-runner scatter-workflow.cwl scatter-job.yml
[workflow scatter-workflow.cwl] start
[step echo] start
[job echo] /tmp/tmp0hqmg400$ echo \
    'Hello world!'
Hello world!
[job echo] completed success
[step echo] start
[job echo_2] /tmp/tmpu65_m1zw$ echo \
    'Hola mundo!'
Hola mundo!
[job echo_2] completed success
[step echo] start
[job echo_3] /tmp/tmp5cs7a2wh$ echo \
    'Bonjour le monde!'
Bonjour le monde!
[job echo_3] completed success
[step echo] start
[job echo_4] /tmp/tmp301wo7p8$ echo \
    'Hallo welt!'
Hallo welt!
[job echo_4] completed success
[step echo] completed success
[workflow scatter-workflow.cwl] completed success
{}
Final process status is success
~~~ 
{: .output}

You can see that the workflow calls echo multiple times on each element of our 
`message_array`. Ok, so how about if we want to scatter over two steps in a workflow?

Let's perform a simple echo like above, but capturing `stdout` by adding the following
lines instead of `outputs: []`

***echo-mod.cwl***

~~~
outputs:
  echo_out:
    type: stdout
~~~
{: .source}

And add a second step that uses `wc` to count the characters in each file. See the tool below:

***wc.cwl***

~~~
{% include cwl/wc.cwl %}
~~~
{: .source}

Now, how do we incorporate scatter? Remember the scatter field is under each step:

***scatter-two-steps.cwl***

~~~
{% include cwl/scatter-two-steps.cwl %}
~~~
{: .source}

Here we have placed the scatter field under each step. This is fine for this example since
it runs quickly, but if you're running many samples for a more complex workflow, you may
wish to consider an alternative. Here we are running scatter on each step independently, but
since the second step is not dependent on the first step completing all languages, we aren't
using the scatter functionality efficiently. The second step expects an array as input from
the first step, so it will wait until everything in step one is finished before doing anything.
Pretend that `echo Hello World!` takes 1 minute to perform, `wc -c` on the output takes 3 minutes
and that `echo Hallo welt!` takes 5 minutes to perform, and `wc` on that output takes 3 minutes.
Even though `echo Hello World!` could finish in 4 minutes, it will actually finish in 8 minutes
because the first step must wait on `echo Hallo welt!`. You can see how this might not scale
well. 

## Visualizing CWL

### [CWL viewer](https://view.commonwl.org/) 
This tool visualises and lists the details of a CWL workflow with its inputs, outputs and steps and packages the files involved into a downloadable Research Object Bundle (zip file with metadata in a manifest), allowing it to be easily viewed and shared.



### [Rabix](https://rabix.io/) 
Power tools for the Common Workflow Language