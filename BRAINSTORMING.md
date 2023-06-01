Brainstorming
---

### Development

- Increase **reproducibility**: store dependencies' versions (not only source code) or even whole environments (images)

- Introduce a **merge request template/checklist** (https://docs.gitlab.com/ee/user/project/description_templates.html , example: page 11 of https://www.forschungsdaten.org/images/e/e4/04-Glaeser-fdm_dumux.pdf [style guide, pass pipeline tests, check new/changed code is covered by test suite, update documentation (doxygen and latex), update CHANGELOG)
- Introduce a release template/checklist (see "To-do list for new releases" in `ROADMAP.md`)

- Consider Gitlab's "Releases" (essentially named commits) instead of creating branches for releases    
  [JSEP: our current system seems fine to me, especially as we want previous versions with back-ported bug fixes to be easily findable.]


### Comunity

- To better open ALF to outside contributions to the code, facilitate authentication to our Gitlab - it should be possible to use, say, ORCID, Github, Gitlab.com accounts (e.g., https://gitlab.hzdr.de/users/sign_in)

- Offer binary packages (.deb, .rpm, etc.), which could generated through pipelines

- Create FAQ through a self-hosted "question & answer(s)" (Stack Exchange like) site: https://www.question2answer.org/
