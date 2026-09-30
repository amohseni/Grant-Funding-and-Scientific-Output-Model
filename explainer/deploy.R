# Deploy the explainer to shinyapps.io. Run from the repository root after
# rsconnect::setAccountInfo(...) has been done once on this machine.
rsconnect::deployApp(appDir = "explainer", appName = "funding-the-gap",
                     appTitle = "Funding the Gap", forceUpdate = TRUE)
