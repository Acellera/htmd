# (c) 2015-2021 Acellera Ltd http://www.acellera.com
# All Rights Reserved
# Distributed under HTMD Software License Agreement
# No redistribution in whole or part
#
import json
import os
import platform
import requests
from subprocess import call
import sys


def show_news():
    try:
        res = requests.get("https://www.htmd.org/news/content", timeout=10)
        print(res.text)
    except:
        pass


def check_approval(product, reg_file):
    reg_data = {}
    try:
        with open(os.path.join(reg_file), "r") as f:
            reg_data = json.load(f)
    except:
        # Couldn't read registration file
        return False

    if not reg_data["ret"]:
        # There's no registration code in the file
        return False

    # Send the registration data
    data = {"code": reg_data["ret"], "product": product}
    url = "https://www.acellera.com/licensing/htmd/check.php"
    res = requests.post(url, data, timeout=10)

    # Check the response
    if res.status_code == 200:
        status = res.json()

        if "approved" in status:
            return True

        if "pending" in status:
            return True

    print('Registration is not approved: %s' % res.text)
    return False


def htmd_registration(product="htmd"):

    reg_file = os.path.join(os.path.expanduser("~"), ".htmd", ".registered-htmd", "registration")

    if not (os.path.isfile(reg_file) and check_approval(product, reg_file)):
        accept_licence(product=product)
        htmd_register(product=product)


def accept_licence(product=None):
    show_licence(product=product)

    ret = input("Type 'yes' to accept the license [no] : ").strip()
    if ret != "yes":
        sys.exit(1)


def show_licence(product):

    # Find the licence file
    path = os.path.join(os.path.dirname(__file__), '..', product, "LICENCE.txt")
    if not os.path.exists(path):
        print("Cannot find the licence file!")
        sys.exit(1)

    # Show the licence file
    if sys.stdout.isatty() and sys.stdin:
        if platform.system() == "Windows":
            call(["c:\\Windows\\System32\\more.com", path])
        else:
            call(["more", path])
    else:
        with open(path) as fh:
            print(fh.read())

def ask(prompt):

    value = ""
    while value == "":
        value = input(prompt).strip()

    return value

def htmd_register(product="htmd"):

    print(
"""
  Welcome to the HTMD registration!

  We would like to know about you to keep in touch.
  Please provide your full name, institutional email,
  institution name, city, and country.
""")

    # Ask a user for data
    data = {}
    data["name"]        = ask("  Full name           : ")
    data["email"]       = ask("  Institutional email : ")
    data["institution"] = ask("  Institution name    : ")
    data["city"]        = ask("  City                : ")
    data["country"]     = ask("  Country             : ")
    data["product"]     = product

    # Send data to the registration server
    url = "https://www.acellera.com/licensing/htmd/register.php"
    res = requests.post(url, data=data, timeout=10)

    # Check the response
    if res.status_code == 200:
        prefix = os.path.join(os.path.expanduser("~"), ".htmd", ".registered-" + product)
        os.makedirs(prefix, exist_ok=True)
        reg_file = os.path.join(prefix, "registration")
        with open(reg_file, "w") as fh:
            fh.write(res.text)
        print("\n  Registration completed!\n")
    else:
        print("\n  Registration failed: %s\n" % res.text)


if __name__ == "__main__":
    pass
