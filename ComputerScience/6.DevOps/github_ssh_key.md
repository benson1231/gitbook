# GitHub SSH Key Setup Guide

This guide provides a clear, step-by-step workflow to generate, register, and use SSH keys with GitHub for secure passwordless authentication.

---

## 🔐 1. Check Existing SSH Keys

Before creating a new key, check whether your system already has one:

```bash
ls ~/.ssh
```

If you see files like:

* `id_rsa` / `id_rsa.pub`
* `id_ed25519` / `id_ed25519.pub`

You already have SSH keys.

---

## 🗝 2. Generate a New SSH Key

GitHub recommends the modern **ed25519** algorithm:

```bash
ssh-keygen -t ed25519 -C "your_email@example.com"
```

Press Enter to accept defaults.

Your keys will be saved to:

* `~/.ssh/id_ed25519`
* `~/.ssh/id_ed25519.pub`

---

## 🔧 3. Start the SSH Agent

Make sure the agent is running:

```bash
eval "$(ssh-agent -s)"
```

Add your private key:

```bash
ssh-add ~/.ssh/id_ed25519
```

---

## 📋 4. Add Public Key to GitHub

Copy the contents of your public key:

```bash
cat ~/.ssh/id_ed25519.pub
```

Go to:
**GitHub → Settings → SSH and GPG Keys → New SSH Key**

Paste the key content.

---

## 🛰 5. Verify SSH Authentication

Test connection:

```bash
ssh -T git@github.com
```

Success message:

```
Hi username! You've successfully authenticated...
```

---

## 🔄 6. Convert Repository Remote URL to SSH

If your repo is using HTTPS, change it:

```bash
git remote set-url origin git@github.com:USERNAME/REPOSITORY.git
```

Check result:

```bash
git remote -v
```

---

## 🚀 7. Push Using SSH

Now pushing won't require a username or password:

```bash
git push
```

---

## 🧰 Optional: SSH Config for Multiple Keys

Create or edit config file:

```bash
nano ~/.ssh/config
```

Example setup:

```
Host github.com
  HostName github.com
  User git
  IdentityFile ~/.ssh/id_ed25519
  IdentitiesOnly yes
```

---

## ✅ Summary

With SSH keys set up, you can:

* Push & pull without typing username/password
* Securely authenticate to GitHub
* Use multiple repos and accounts smoothly

Your Git workflow becomes faster and more secure.

---

If you'd like, I can also add:

* Windows + WSL2 SSH setup
* Multi-account GitHub SSH setup
* Auto script to convert all local repos to SSH URLs.
