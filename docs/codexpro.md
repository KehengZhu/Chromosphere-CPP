# CodexPro

This project has a saved CodexPro workspace profile for:

```text
/Users/zkeheng/SWMFSoftware/Chromosphere2026
```

The profile uses a dedicated Cloudflare named tunnel:

```text
https://codex-chromosphere2026.kehengphysics.site
```

The saved profile lives outside the repository under `~/.codexpro/profiles/` and contains the private CodexPro token. Do not copy that token into this repository.

## Daily Start

From the project root:

```bash
codexpro start
```

or, equivalently:

```bash
codexpro stable \
  --root /Users/zkeheng/SWMFSoftware/Chromosphere2026 \
  --hostname codex-chromosphere2026.kehengphysics.site \
  --tunnel-name codex-chromosphere2026
```

Leave that terminal running while ChatGPT is connected.

## ChatGPT Plugin Connector

Use ChatGPT Developer Mode, then create or update the developer-mode plugin for the CodexPro MCP server.

In ChatGPT:

```text
Settings -> Security and login -> Developer mode: on
Settings -> Plugins -> + button
```

Use:

```text
Connection: Server URL
Authentication: None / No Authentication
```

The server URL has this shape:

```text
https://codex-chromosphere2026.kehengphysics.site/mcp?codexpro_token=...
```

CodexPro prints and copies the full URL, including the private token, when it starts.

Replace the old ngrok Server URL in the plugin with this newly printed URL. The token is private; do not paste it into project files or chat messages.

## Current Profile Settings

```text
Tunnel: cloudflare-named
Tunnel name: codex-chromosphere2026
Hostname: codex-chromosphere2026.kehengphysics.site
Local port: 8793
Mode: agent
Bash: safe
Write: workspace
Tool mode: standard
.gitignore: on
Codex sessions: off
```

## Useful Commands

Check the saved profile:

```bash
codexpro settings show
```

Check runtime readiness:

```bash
codexpro doctor
```

Re-save the stable Cloudflare profile:

```bash
codexpro settings set \
  --root /Users/zkeheng/SWMFSoftware/Chromosphere2026 \
  --port 8793 \
  --tunnel cloudflare-named \
  --hostname codex-chromosphere2026.kehengphysics.site \
  --tunnel-name codex-chromosphere2026 \
  --mode agent \
  --bash safe \
  --write workspace \
  --tool-mode standard \
  --honor-gitignore
```

Remove this workspace's saved CodexPro profile:

```bash
codexpro settings delete --yes
```
