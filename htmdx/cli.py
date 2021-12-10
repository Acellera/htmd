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


def htmd_show_news():
    try:
        res = requests.get("https://www.htmd.org/news/content", timeout=3)
        print(res.text)
    except:
        pass


def _get_reg_file(product):

    dir = os.path.join(os.path.expanduser("~"), "." + product, ".registered-" + product)
    os.makedirs(dir, exist_ok=True)
    file = os.path.join(dir, "registration")

    return file


def _check_registration(product):

    # Get a registation file
    reg_file = _get_reg_file(product)

    if not os.path.exists(reg_file):
        print("Registation file does not exist")
        return False

    data = {"product": product}
    try:
        with open(reg_file) as fh:
            data.update(json.load(fh))
    except:
        print("Cannot read the registration file")
        return False

    # Note: remove after switching to the new backend
    if "code" not in data and "ret" in data:
        data["code"] = data["ret"]

    # Send the registration data
    url = "https://www.acellera.com/registration/check"
    res = requests.post(url, data, timeout=10)

    # Check the response
    if res.status_code == 200:
        status = res.json()
        if "approved" in status:
            return True

        # Note: remove after switching to the new backend
        if "pending" in status:
            return True

    print(f"Registration is not approved: {res.text}")
    return False


def _accept_licence(product):

    ret = input("Type 'yes' to accept the license [no] : ").strip()
    if ret != "yes":
        sys.exit(1)


def _show_licence(product):

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


def htmd_registration(product="htmd"):

    if not _check_registration(product):
        _show_licence(product)
        _accept_licence(product)
        htmd_register(product)


def _ask(prompt):

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
    data["name"]        = _ask("  Full name           : ")
    data["email"]       = _ask("  Institutional email : ")
    data["institution"] = _ask("  Institution name    : ")
    data["city"]        = _ask("  City                : ")
    data["country"]     = _ask("  Country             : ")
    data["product"]     = product

    # Send data to the registration server
    url = "https://www.acellera.com/registration/register"
    res = requests.post(url, data=data, timeout=10)

    # Check the response
    if res.status_code == 200:
        reg_file = _get_reg_file(product)
        with open(reg_file, "w") as fh:
            fh.write(res.text)
        print("\n  Registration completed!\n")
    else:
        print(f"\n  Registration failed: {res.text}\n")


if __name__ == "__main__":
    pass
