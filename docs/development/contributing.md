# Development conventions

Follow the root CONTRIBUTING.md and preserve licenses. Use snake_case Python modules/functions, PascalCase classes and UPPER_SNAKE_CASE constants. New public modules need semantic names, an existing registry mapping and an explicit protocol. Historical aliases are not a naming template.

Audit before modifying scientific code. Keep experimental, supported, deprecated and archived paths distinct. Do not replace an engine default or remove a check to pass tests. Keep user changes, historical results and current running tasks isolated. Never force-push, auto-merge, upgrade scientific engines or publish a stable release as part of a governance-only refactor.

For a new algorithm, add a backend under the existing functional module and a separate scientific protocol when warranted. Supply dependency and environment checks, exact input/output identity, compatibility mapping and tests. Source presence, installation, execution validation and biological evidence remain separate.
