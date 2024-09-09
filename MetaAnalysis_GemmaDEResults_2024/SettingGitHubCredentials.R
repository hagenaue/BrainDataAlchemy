#Setting GitHub Credentials
#Megan Hagenauer, July 1, 2024

#Checking on my GitHub settings:
library("usethis")
git_sitrep()

#Example Output:
# ── Git global (user) 
# • Name: 'Megan Hastings Hagenauer'
# • Email: 'hagenaue@umich.edu'
# • Global (user-level) gitignore file: <unset>
#   • Vaccinated: FALSE
# ℹ See `?git_vaccinate` to learn more
# ℹ Defaulting to 'https' Git protocol
# • Default Git protocol: 'https'
# • Default initial branch name: <unset>
#   
#   ── GitHub user 
# • Default GitHub host: 'https://github.com'
# • Personal access token for 'https://github.com': <unset>
#   • To create a personal access token, call `create_github_token()`
# • To store a token for current and future use, call `gitcreds::gitcreds_set()`
# ℹ Read more in the 'Managing Git(Hub) Credentials' article:
#   https://usethis.r-lib.org/articles/articles/git-credentials.html
# 
# ── Active usethis project: '/Users/hagenaue/Documents/Example2024_BrainDataAlchemy' ──
# 
# ── Git local (project) 
# • Name: 'Megan Hastings Hagenauer'
# • Email: 'hagenaue@umich.edu'
# • Default branch: 'main'
# • Current local branch -> remote tracking branch:
#   'main' -> 'origin/main'
# 
# ── GitHub project 
# • Type = 'maybe_ours_or_theirs'
# • Host = 'https://github.com'
# • Config supports a pull request = NA
# • origin = 'hagenaue/Example2024_BrainDataAlchemy'
# • upstream = <not configured>
#   • Desc = 'origin' is a GitHub repo and 'upstream' is either not configured or is not a GitHub repo.
# 
# We may be offline or you may need to configure a GitHub personal access
# token. `gh_token_help()` can help with that.
# 
# Read more about what this GitHub remote configurations means at:
#   'https://happygitwithr.com/common-remote-setups.html'

create_github_token()

gitcreds::gitcreds_set()
# -> Adding new credentials...
# -> Removing credentials from cache...
# -> Done.
