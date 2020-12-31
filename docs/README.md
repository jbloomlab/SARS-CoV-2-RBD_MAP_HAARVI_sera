# How to set up a landing page for `dms-view` links
SKH

Here are steps to set up a simple webpage using `GitHub pages`.

The URL is `<username>.github.io/<repo name>` so for this specific project the URL is [https://jbloomlab.github.io/SARS-CoV-2-RBD_MAP_HAARVI_sera/](https://jbloomlab.github.io/SARS-CoV-2-RBD_MAP_HAARVI_sera/).

## 1. Create a docs directory in your repo

## 2. Create the neccessary files for `GitHub pages`.

There are two files you need to create the website.

[`./index.md`](./index.md): A markdown file where you will create the content for the website.

[`./_config.yml`](./_config.yml): A yaml file to set some website attributes. Make sure you change the information to match the information of your new project. *importantly*, the baseurl field should be should be "/repo_name" and the url field should be an empty string, "".

## 3. Turn on `GitHub pages`

You should specify that the website is served through the `docs` repo from the `master` branch.

## 4. Edit your content

### adding text

To add text to the landing page, edit the [`index.md`](index.md) file.
You can write the text using markdown syntax, including links to `dms-view` pages.
Push the changes to the `master` branch and they will show up automatically on the website (be patient, sometimes it takes a few minutes for the webpage to update).

## changing the title

If you would like to change the title of the webage, edit line 21 `title: ` of the [`_config.yml`](_config.yml) file.
