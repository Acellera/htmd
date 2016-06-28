# HTMD Versions Explained

HTMD versions are of the following format:

```
<big release>.<major release (stable and develment tagging)>.<minor release (bug fixes)>
```

### Big release

Only changed at the developers' discretion. This may only happen when there are ground-breaking changes to the software.

### Major release - Definition of _stable_ and _development_ versions

Major releases are periodic (every 3 or 6 months, check the [milestones](https://github.com/Acellera/htmd/milestones) for when the next one comes out). When a release is done, actually two versions come out: a __stable__ one and a __development__ one.

When the second/middle number (major release) is even (0,2,...,2n), it means that this major release is a __stable__ release. Examples: 1.0.0; 1.2.0; 1.24.0.
In this context, the third/right number will be 0 (zero) when the major release comes out and will increase until a new major release is done. 

For each stable release, there is a corresponding __development__ release, whose second/middle number is odd (1,3,...,2n+1) and always one more than the stable. Examples: 1.1.0; 1.3.0; 1.25.0.
The third/right number is also a 0 (zero) when the major release comes out and will increase until a new major release is done.

### Minor release

Minor releases correspond to the third/right number and are, in general, periodic (every two weeks). In the case of being a _stable_ version, they correspond to critical bug-fixes and may be released as soon as the fix is made. 

This number is always 0 (zero) at major releases and the corresponding _stable_ and _development_ versions (example: 1.0.0 and 1.1.0 correspond) are "in sync" when this happens.

When this number is different from 0 (zero), nothing should be assumed in the relationship between corresponding _stable_ and _development_ versions.

### GitHub Tagging

Stable releases are tagged on a specific branch created for each of these (named `rel-<big release>.<major release>.x`). Development releases are tagged on the master branch.

Useful tag listing commands:
```
git tag -n
git describe --tags
```

# How to release a new HTMD version?

These are examples of how to release HTMD versions.

DO NOT ONLY do `git push --tags`, as it will push the tag but not the commit!

### Big and/or major releases

Imagine one wants to do a big/major release (in this case, let's assume it's the first major release of big release 1).

1. Make sure you are working on `Acellera/htmd:master` (https://github.com/Acellera/htmd.git):

   ```
   git remote -v
   git checkout master
   ```

1. Make sure the `master` branch is up-to-date:

   ```
   git fetch
   git pull
   ```

1. On `master`, create the new stable branch, check it out, and tag it:

   ```
   git branch rel-1.0.x
   git checkout rel-1.0.x
   git tag -a 1.0.0 -m "new stable release"
   ```

1. Push the new branch and tag to the remote (`origin`):

   `git push --tags origin rel-1.0.x`

   This will trigger two Travis builds: one due to the branch and another due to the tag. A conda release will be made.
1. Check out `master`, tag the development version, and push the tag:

   ```
   git checkout master
   git tag -a 1.1.0 -m "new development release"
   git push --tags
   ```

   This will trigger another Travis build. A conda release will be made.

These two tags will point to the same commit (in sync).

### Minor releases - stable

Imagine one wants to do a minor release (bug-fix) on release 1.0.0.

1. Make sure you are working on `Acellera/htmd:rel-1.0.x` (https://github.com/Acellera/htmd.git):

   ```
   git remote -v
   git fetch
   git checkout rel-1.0.x
   ```

1. Make sure the `rel-1.0.x` branch is up-to-date:

   ```
   git fetch
   git pull
   ```

1. Do the fix, add the files, and commit it.
1. Tag the new bug-fix:

   `git tag -a 1.0.1 -m "new bug-fix"`

1. Push the fix and the tag to the remote (`origin`):

   `git push --tags origin rel-1.0.x`

   This will trigger a new Travis build. A conda release will be made.

In alternative, push the commit but do not tag it, if it is not critical to release right now.

### Minor releases - development

Imagine one wants to do a minor release on release 1.1.0.

1. Make sure you are working on `Acellera/htmd:master` (https://github.com/Acellera/htmd.git):

   ```
   git remote -v
   git checkout master
   ```

1. Make sure the `master` branch is up-to-date:

   ```
   git fetch
   git pull
   ```

1. Tag the new minor release:

   `git tag -a 1.1.1 -m "new minor release"`

1. Push the tag to the remote (`origin`):

   `git push --tags origin master`

   This will trigger a new Travis build. A conda release will be made.
