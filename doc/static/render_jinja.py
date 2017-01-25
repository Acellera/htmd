from jinja2 import Environment, FileSystemLoader
import os
from glob import glob
import re
import sys

ind = sys.argv[1]
outd = sys.argv[2]

def render_html(indir, outdir):
    indir = os.path.abspath(indir)
    outdir = os.path.abspath(outdir)

    templatedir = os.path.join(indir, 'templates')

    env = Environment(loader=FileSystemLoader(templatedir))
    for t in env.list_templates():
        if t == 'common.html':
            continue
        temp = env.get_template(t)
        html = temp.render({'name': t})
        with open(os.path.join(outdir, t), 'w') as f:
            f.write(html)

render_html(ind, outd)
    

