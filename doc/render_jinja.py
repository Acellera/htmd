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
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    env = Environment(loader=FileSystemLoader(indir))
    for t in env.list_templates():
        if t == 'common.html' or t == 'navbar.html':
            continue
        temp = env.get_template(t)
        html = temp.render()
        with open(os.path.join(outdir, t), 'w') as f:
            f.write(html)

render_html(ind, outd)
    

