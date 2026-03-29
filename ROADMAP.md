# Roadmap


Merges
---

* Development branches to be merged to master only when their 
  - **documentation** reflects all changes
  - **CHANGELOG.md** mentions all breaking changes and critical bug fixes


Release cycle
---

* Half-yearly releases: **April and October**
* Main release (new documentation paper): ~**every three to four years**

Latest main release documentation: ALF 2.4
Former: ALF 2.0, published on 2022 (2020 on ArXiV).


 To-do list for new releases:
 ---
  - check the pipeline tests pass
  - check the main documentation (paper) content is up-to-date (i.e., all code changes are reflected in the text)
  - check CHANGELOG is up-to-date
  - create `ALF-#.#` release branch from master, then
    - update main documentation (paper — set LaTeX variables such as `\ALFver` appropriately)
    - update README file:
      - add version number to header of  (ALF -> ALF #.#)
      - remove warning: "-> This is the development version, the latest..."
  - update **master**'s README warning: "-> This is the development version of ALF, the latest stable version is ALF #.#"
  - create corresponding pyALF branch (or add remark to its README it's compatible with the new version)
  - create corresponding Tutorial branch (or add remark to its README it's compatible with the new version)
  - update ALF website (links, and also 'news') and advertise in other channels (e.g., Discord announcements)


Next Goals
---

(To be updated)

- Expand test pipelines to include simulations
- Improved overall usability:
  - Further develop pyALF
  - **Allow models to be input instead of coded**
- Time-dependent Hamiltonians
- Hybrid Monte Carlo with exact forces
- Interaction expansion (CT-INT)

**Check ideas in `BRAINSTORMING.md`**.


ALF 2.0
---
(_former "ALF 1.2"_)

- Parallel tempering
- Global updates / global tau
- Projective approaches
- Continuous fields
- Langevin
- Maxent
- More models
- Trotter symmetric
- Predefined structures.
The last two still need some work from Fakher, and Langevin has to be merged in, but then the code will be essentially ready for release and we only have to worry about the documentation (paper).


