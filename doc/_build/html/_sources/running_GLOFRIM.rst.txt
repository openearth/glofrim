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

.. sourcecode:: ipython

   In [1]: x = 2

   In [2]: x**3
   Out[2]: 8

When you reload the page by refreshing your browser pointing to
:file:`_build/html/index.html`, you should see a link to the
"Getting Started" docs, and in there this page with the screenshot.
`Voila!`

We can also use the image directive in :file:`running_GLOFRIM.rst` to include to the screenshot above
with::

  .. image::
     _static/GermanDischarge.png



   
