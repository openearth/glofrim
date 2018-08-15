.. _running_GLOFRIM:


***************
Running GLOFRIM
***************

Installing your doc directory
=============================

GLOFRIM in da house whoop whoop what's up and hello World!

You may already have sphinx `sphinx <http://sphinx.pocoo.org/>`_
installed -- you can check by doing::

  python glofrim_runner.py run /path/to/glofrim.ini

Or otherwise try something else.

You can also add an image:

.. image:: _static/GermanDischarge.png

.. _fetching-the-data:

Fetching the data
-----------------

Now we will start to customize out docs.  Grab a couple of files from
the `web site <https://github.com/matplotlib/GLOFRIM>`_
or git.  You will need :file:`getting_started.rst` and
:file:`_static/basic_screenshot.png`.  All of the files live in the
"completed" version of this tutorial, but since this is a tutorial,
we'll just grab them one at a time, so you can learn what needs to be
changed where.  Since we have more files to come, I'm going to grab
the whole git directory and just copy the files I need over for now.
First, I'll cd up back into the directory containing my project, check
out the "finished" product from git, and then copy in just the files I
need into my :file:`GLOFRIM` directory::

  home:~/tmp/GLOFRIM> pwd
  /Users/jdhunter/tmp/GLOFRIM
  home:~/tmp/GLOFRIM> cd ..
  home:~/tmp> git clone https://github.com/matplotlib/GLOFRIM.git tutorial
  Cloning into 'tutorial'...
  remote: Counting objects: 87, done.
  remote: Compressing objects: 100% (43/43), done.
  remote: Total 87 (delta 45), reused 83 (delta 41)
  Unpacking objects: 100% (87/87), done.
  Checking connectivity... done
  home:~/tmp> cp tutorial/getting_started.rst GLOFRIM/
  home:~/tmp> cp tutorial/_static/basic_screenshot.png GLOFRIM/_static/

The last step is to modify :file:`index.rst` to include the
:file:`getting_started.rst` file (be careful with the indentation, the
"g" in "getting_started" should line up with the ':' in ``:maxdepth``::

  Contents:

  .. toctree::
     :maxdepth: 2

     getting_started.rst

and then rebuild the docs::

  cd GLOFRIM
  make html

.. sourcecode:: ipython

    In [69]: lines = plot([1,2,3])

    In [70]: setp(lines)
      alpha: float
      animated: [True | False]
      antialiased or aa: [True | False]
      ...snip

When you reload the page by refreshing your browser pointing to
:file:`_build/html/index.html`, you should see a link to the
"Getting Started" docs, and in there this page with the screenshot.
`Voila!`

We can also use the image directive in :file:`index.rst` to include to the screenshot above
with::

  .. image::
     _static/GermanDischarge.png
