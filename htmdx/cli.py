# (c) 2015-2018 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import os
import sys
from subprocess import call
import inspect
import htmdx
import json
import requests
import platform

has_connection = True


def main_do_nothing():
    import htmd


def htmd_do_register():
    do_register('htmd')

def show_news():
    try:
        r = requests.get("https://www.htmd.org/news/content", timeout=(3.05, 3.05))
        print(r.content.decode("ascii"))
    except:
        pass


def check_approval(product, reg_file):
    # TODO: this may come back in other form
    # from htmdx.license import licenseEntitlements
    # jj = licenseEntitlements()
    # if "HTMD" in jj: return True
    # if "htmd" in jj: return True

    j = {}
    try:
        with open(os.path.join(reg_file), "r") as f:
            j = json.load(f)
    except:
        # Couldn't read registration file
        return False

    if not j['ret']:
        # There's no registration code in the file
        return False

    # Do online check
    try:
        payload = {
            'code': j['ret'],
            'product': product
        }
        r = requests.post("https://www.acellera.com/licensing/htmd/check.php", params=payload, timeout=(3.05, 3.05))
        ret = json.loads(r.content.decode("ascii"))
        if 'approved' in ret:
            return True
        if 'pending' in ret:
            return True
        #    except requests.exceptions.Timeout as e:
        # dont' trigger registration if connection timed out
        #      return True
    except:
        return True

    return False


def check_registration(product=None):
    if not product:
        product = "NA"

    reg_file = os.path.join(os.path.expanduser('~'), '.htmd', '.registered-htmd', 'registration')

    if os.path.isfile(reg_file):
        if not check_approval(product, reg_file):
            accept_license(product=product)
            do_register(product=product)
    else:
        accept_license(product=product)
        do_register(product=product)


def accept_license(product=None):
    show_license(product=product)
    print("Type 'yes' to accept the license [no] : ", end="")
    sys.stdout.flush()
    ret = input().strip()
    if ret != "yes":
        os._exit(1)


def show_license(product=None):
    path = os.path.abspath(os.path.dirname(os.path.dirname(inspect.getfile(accept_license))))
    if (product):
        path = os.path.join(path, product)

    path = os.path.join(path, "LICENCE.txt")
    if not os.path.exists(path):
        print("Could not find license file. Aborting")
        os._exit(1)

    print(sys.stdin)
    if sys.stdout.isatty() and sys.stdin:
        if platform.system() == "Windows":
           call(["c:\\Windows\\System32\\more.com", path])
        else:
           call(["more", path])
    else:
        fh = open(path, "r")
        print(fh.read())
        fh.close()


def do_register(product=None):
    name = ""
    institution = ""
    email = ""
    city = ""
    country = ""

    print("\n Welcome to HTMD. We'd like to know who you are so that we can keep in touch!")
    print(" Please enter your name, affiliation and contact email address.\n")

    while email == "" or not ("@" in email):
        print(" Email       : ", end="")
        sys.stdout.flush()
        email = input().strip()

    while name == "":
        print(" Full name   : ", end="")
        sys.stdout.flush()
        name = input().strip()
    while institution == "":
        print(" Institution : ", end="")
        sys.stdout.flush()
        institution = input().strip()

    while city == "":
        print(" City        : ", end="")
        sys.stdout.flush()
        city = input().strip()

    while country == "":
        print(" Country     : ", end="")
        sys.stdout.flush()
        country = input().strip()

    import requests
    payload = {
        'name': name,
        'institution': institution,
        'email': email,
        'product': product,
        'city': city,
        'country': country
    }
    try:
        r = requests.post("https://www.acellera.com/licensing/htmd/register.php", params=payload)
        ret = json.dumps(r.content.decode("ascii"))

        prefix = os.path.join(os.path.expanduser('~'), '.htmd')
        if not os.path.exists(prefix):
            os.mkdir(prefix)
        prefix = os.path.join(prefix, '.registered-' + product)
        if not os.path.exists(prefix):
            os.mkdir(prefix)
        regfile = os.path.join(prefix, "registration")
        fh = open(regfile, "w")
        fh.write(r.content.decode("ascii"))
        fh.close()
        print("")
    except:
        print("\nCould not register!\n")


def check_ipython_profile():
    prefix = os.path.join(os.path.expanduser('~'), '.ipython')
    if not os.path.exists(prefix):
        os.mkdir(prefix)
    prefix = os.path.join(prefix, 'profile_default')
    if not os.path.exists(prefix):
        os.mkdir(prefix)
    prefix = os.path.join(prefix, 'startup')
    if not os.path.exists(prefix):
        os.mkdir(prefix)

    prefix = os.path.join(prefix, "00-htmd.py")
    if not os.path.exists(prefix):
        f = open(prefix, "w")
        f.write("from htmd.ui import *\n")
        f.close()


def main_activate():
    import sys
    import requests
    if len(sys.argv) != 2:
      print("\n Acellera License Installation\n\n Syntax: activate [activation token]\n\nPlease contact support@acellera.com for licensing\n\n" )
      sys.exit(0)
    
    token=sys.argv[1]
    r = requests.get( "https://licensing.acellera.com/issue/?token=" + token )
    if ( r.status_code != 200 ):
       print("\n Error: " + r.text + "(" + str(r.status_code) + ")\n\n" )
    else:
       prefix = os.path.join(os.path.expanduser('~'), ".acellera" );
       try: 
         os.mkdir( prefix )
       except:
         pass
       f = open(os.path.join( prefix, "license.dat"), "a" )   
       try:      
         print(r.content.decode("ascii"), file=f )
       except:
         print(r.content, file=f )
       f.close()
       print( "\n License installed in ~/.acellera/license.dat\n\n" );
       try:
         f = open(os.path.join( "/opt/acellera/license.dat"), "a" )   
            
         print(r.content.decode("ascii"), file=f )
         f.close()
         print( "\n License installed in /opt/acellera/license.dat\n\n" );
       except:
         pass


def main_htmd():
    check_registration(product="htmd")

    #    check_ipython_profile()
    sys.argv.append('--no-banner')
    sys.argv.append('--no-confirm-exit')
    if ('DISPLAY' in os.environ) and os.environ['DISPLAY']:
        import qtconsole.qtconsoleapp
        qtconsole.qtconsoleapp.main()
    else:
        import jupyter_console.app
        jupyter_console.app.main()


def main_htmd_notebook():
    import sys
    check_registration(product="htmd")

    #    check_ipython_profile()

    from notebook.notebookapp import main
    sys.exit(main())


if __name__ == '__main__':
    main_activate()
