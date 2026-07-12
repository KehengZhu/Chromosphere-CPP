# CodexPro

This project has a saved CodexPro workspace profile for:

```text
/Users/zkeheng/SWMFSoftware/Chromosphere2026
```

The profile uses the ngrok stable dev domain:

```text
unlatch-doorpost-very.ngrok-free.dev
```

The saved profile lives outside the repository under `~/.codexpro/profiles/` and contains the private CodexPro token. Do not copy that token into this repository.

## Daily Start

From the project root:

```bash
codexpro start
```

or, equivalently:

```bash
codexpro ngrok \
  --root /Users/zkeheng/SWMFSoftware/Chromosphere2026 \
  --hostname unlatch-doorpost-very.ngrok-free.dev \
  --bash safe
```

Leave that terminal running while ChatGPT is connected.

## ChatGPT Plugin Connector

Use ChatGPT Developer Mode, then create a developer-mode plugin/app entry for the CodexPro MCP server.

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
https://unlatch-doorpost-very.ngrok-free.dev/mcp?codexpro_token=...
```

CodexPro prints and copies the full URL, including the private token, when it starts.

If ChatGPT uses the word "Plugins" instead of "Apps", that is expected. The developer-mode plugin entry is still the place where the CodexPro MCP Server URL goes.

## Current Profile Settings

```text
Tunnel: ngrok
Hostname: unlatch-doorpost-very.ngrok-free.dev
Local port: 8787
Mode: agent
Bash: safe
Write: workspace
Tool mode: standard
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

Change the saved ngrok hostname:

```bash
codexpro settings set --tunnel ngrok --hostname NEW-HOSTNAME.ngrok-free.dev
```

Remove this workspace's saved CodexPro profile:

```bash
codexpro settings delete --yes
```
