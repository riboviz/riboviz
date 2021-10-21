# Writing and updating documentation

* [Updating workflow images](#updating-workflow-images)

---

## Updating workflow images

Workflow images in `docs/images/` are written in the open source [dot](https://graphviz.org/doc/info/lang.html) language from [GraphViz](https://www.graphviz.org/). For an overview, see [Drawing graphs with dot](https://www.graphviz.org/pdf/dotguide.pdf).

If you update the `.dot` files, you should also update the corresponding `.svg` images. GraphViz includes a command-line tool, `dot`, for converting dot files into image in various formats.

To convert a `.dot` file to an `.svg` file:

```console
$ dot -Tsvg workflow.dot > workflow.svg
```

Alternatively, to convert all `.dot` files in `docs/images/`, use the `Makefile`:

```console
$ make clean svgs
```

Alternatively, use an [Editor supporting live preview of GraphViz images](./install.md#editor-supporting-live-preview-of-graphviz-images-optional) that also allows these to be exported as `.svg` images.
