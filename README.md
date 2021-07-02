# PyDBHelper
minimal coverage, third normal form, and much more

It works in two modes: the first one is for schemas, the second one is for physical data organization.

Note: the code is a mess, but I wasn't looking for performance and perfection.

**For schemas:**
- tells if a scheme is in 3NF
- calculates minimal coverage
- find primary keys
- given a set of keys, it tells if they're key or not
- closure of a set of attributes

**For physical organization data:**
- solve isam exercises:
  1. Size of Bucket Directory and hash file
  2. Search cost
  3. Number of buckets to have a search cost less equal than <number provided>

- solve hash exercises

**WIP functionalities:**
- check if decompositions in ro preserve FDs --> currently it states the FDs you manually need to check
- check if R, given decompositions in ro, has a lossless join. --> currently it prints the big table at its initial state
- solve btree exercises

**Help is appreciated!**

## Requirements
Python3 installed on your machine

## Instructions
Download or clone this repo inside a folder.

Run `python DBHelper.py` for help (shown below)

```
Syntax: DBHelper.py <option> <parameter:json_file_containing_a_scheme | all>
<option>:
        -h --help:
                displays this message.
        -v --verbose:
                shows `what's going on` messages
        -x --extra-verbose:
                shows anything that this program could print
        -n --third-nf:
                tells if a scheme is in 3NF or not
        -k --find-keys:
                find and shows keys of a scheme
        -t --test-keys:
                test keys defined in the scheme, if they are keys or not
        -c --closure:
                closure of a given X attributes
        -m --minimal:
                find a minimal coverage of a scheme
        [WIP] -p --preserve:
                check if decompositions in ro preserve FDs --> currently it states the FDs you manually need to check
        [WIP] -l --lossless:
                check if R, given decompositions in ro, has a lossless join. --> currently it prints the big table at its initial state.
        --json:
                custom json folder (default: .)

PHYSICAL DATA ORGANIZATION
        --physical-isam:
                solve a json containing isam exercise
        --physical-hash:
                solve a json containing hash exercise
        [WIP]--physical-btree:
                solve a json containing btree exercise
        --phyjson:
                custom json folder (default: phyjson)
NOTE: WIP functions are WIP for a reason: you're welcome to give a hand where needed!

```

# JSON Structures
There are two folders:
- _json_: contains **schemas**
- _phyjson_: contains **physical data specs**

In json folder, specify the following fields (in a named.json file):
- "R": optional but suggested, set of attributes which schema contains
- "F": **mandatory**, functional dependencies
- "X": set of attributes, mandatory if you want to test closure (`--closure` option)
- "ro": set of decompositions, mandatory if you want to use `--lossless` option
- "keys": set of attributes, mandatory if you want to test keys (`--test-keys` option)

Example:
```
{
    "R": "ABCD",

    "F": [
        "AB -> CD",
        "BC -> A",
        "D -> AC"
    ],
    "X": "AB",
    "ro": ["ACD, ABE, CDG"],
    "keys": "AB, BC, BD"
}
```

In the `phyjson` folder, specify the following fields (in a named.json file):
  
**For ISAM**
- "type": it could be ISAM, HASH, or BTREE, it specifies the data type you are going to use
- "NR": total number of records into the file
- "R": single record size
- "K": the size of the key
- "CB": capacity of a block, usually 2048 byte
- "P": pointer size
- "filled_space": how much blocks are filled in percentage (100 means, all) 
- "floor_min": optional field, if set to 1 it will give the minimum number of records in the blocks, by default it gives the maximum

**For HASH**
See above, with these additional fields:
- "no_access": optional, number of accesses, it is <number provided> to calculate the number of buckets to have a search cost less equal than <number provided>
- "B": number of buckets in the data structure
  
Example:
```
{
    "type": "ISAM",
    "NR": 817000,
    "R": 103,
    "K": 31,
    "CB": 2048,
    "P": 5,
    "filled_space": 80,
    "floor_min": 0
}
```
