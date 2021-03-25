# (c) 2015-2021 Acellera Ltd http://www.acellera.com
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
import urllib.request
import urllib.parse
import platform


def htmd_do_register():
    do_register("htmd")


def show_news():
    try:
        response = urllib.request.urlopen(
            "https://www.htmd.org/news/content", timeout=3.05
        )
        print(response.read().decode("ascii"))
    except:
        pass


def check_approval(product, reg_file):
    registration_data = {}
    try:
        with open(os.path.join(reg_file), "r") as f:
            registration_data = json.load(f)
    except:
        # Couldn't read registration file
        return False

    if not registration_data["ret"]:
        # There's no registration code in the file
        return False

    # Do online check
    try:
        payload = {"code": registration_data["ret"], "product": product}
        data = urllib.parse.urlencode(payload).encode("ascii")
        with urllib.request.urlopen(
            "https://www.acellera.com/licensing/htmd/check.php", data, timeout=3.05
        ) as f:
            ret = json.loads(f.read().decode("ascii"))

        if "approved" in ret:
            return True
        if "pending" in ret:
            return True
    except:
        return True

    return False


def check_registration(product=None):
    reg_file = os.path.join(
        os.path.expanduser("~"), ".htmd", ".registered-htmd", "registration"
    )

    if not (os.path.isfile(reg_file) and check_approval(product, reg_file)):
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
    path = os.path.abspath(
        os.path.dirname(os.path.dirname(inspect.getfile(accept_license)))
    )
    if product:
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

    print(
        "\n Welcome to HTMD. We'd like to know who you are so that we can keep in touch!"
    )
    print(" Please enter your name, affiliation and institutional email address.\n")
    print(" Please provide valid information so that it might be verified.\n")

    while email == "" or not ("@" in email) or not ("." in email):
        print(" Institutional Email : ", end="")
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

    payload = {
        "name": name,
        "institution": institution,
        "email": email,
        "product": product,
        "city": city,
        "country": country,
    }
    try:
        data = urllib.parse.urlencode(payload).encode("ascii")
        with urllib.request.urlopen(
            "https://www.acellera.com/licensing/htmd/register.php", data
        ) as f:
            text = f.read().decode("ascii")

        prefix = os.path.join(
            os.path.expanduser("~"), ".htmd", ".registered-" + product
        )
        os.makedirs(prefix, exist_ok=True)
        regfile = os.path.join(prefix, "registration")
        with open(regfile, "w") as fh:
            fh.write(text)
        print("")
    except:
        print("\nCould not register!\n")


if __name__ == "__main__":
    pass
